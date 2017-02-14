from builtins import object
from builtins import str

from bauhaus.resources import getResourcePath
from bauhaus.utils import mkdirp, readFile, writeFile

from collections import OrderedDict, namedtuple
from contextlib import closing, contextmanager
import ninja_syntax as ninja, shutil, os.path as op, jinja2



Rule           = namedtuple("Rule", ("name", "command"))
BuildStatement = namedtuple("BuildStatement", ("outputs", "rule", "inputs", "variables"))

class IncorrectDuplicateRule(Exception): pass

class ContextTracker(object):
    """
    A class that enables tracking a stack of "context" variables
    """
    def __init__(self):
        self._context = []

    def lookup(self, key):
        for (k, v) in self._context:
            if k == key:
                return v
        else:
            raise KeyError(key)

    def __getitem__(self, key):
        return self.lookup(key)

    def _push(self, key, value):
        self._context.insert(0, (key, value))

    def _pop(self):
        (k, v) = self._context[0]
        del self._context[0]
        return (k, v)

    @contextmanager
    def context(self, key, value):
        """
        RAII method for use with "with" statement
        """
        self._push(key, value)
        yield
        self._pop()

    def contextToDict(self):
        return dict(self._context)



BundledResource = namedtuple("BundledResource", ("resourceName", "destPath", "substitutions"))

class PFlow(ContextTracker):

    def __init__(self, logDir=""):
        super(PFlow, self).__init__()
        self._rules = OrderedDict()
        self._buildStmts = []
        self._resourcesToBundle = []
        self.bundleResource("run.sh")
        self._grid = False
        self._nChunks = 8

    # ----- script/resource bundling ---------

    def bundleResource(self, resourceName, destPath=None, substitutions=dict()):
        self._resourcesToBundle.append(BundledResource(resourceName, destPath, substitutions))


    # ----- configuration ------
    # --- TODO: I don't think this belongs in the PFlow class

    def noGrid(self) :
        """Disable the qsub/farm option"""
        self._grid = False

    def grid(self, queueName):
        self._grid = True
        self._sgeQueue = queueName

    @property
    def chunks(self):
        """
        Set the nChunks for workflows that support scatter/gather
        - 0 means "don't do chunks"
        - n>0 means split datasets n ways for scatter/gather

        TODO: I think we ideally want to specify the size of a chunk, not the nChunks
        """
        return self._nChunks

    @chunks.setter
    def chunks(self, nChunks):
        if (not isinstance(nChunks, int)) or (nChunks < 0):
            raise ValueError("Illegal nChunks")
        self._nChunks = nChunks

    # ---- rules, build targets  -----

    def formatInContext(self, s):
        if isinstance(s, (str, bytes)):
            return s.format(**self.contextToDict())
        else:
            return s

    def genRuleOnce(self, name, command):
        # Redundant rules get coalesced; error on specifying rules with
        # same name but different content
        if name in self._rules:
            if self._rules[name] != command:
                raise IncorrectDuplicateRule(name)
        else:
            self._rules[name] = command
        return self._rules[name]

    def genBuildStatement(self, outputs, rule, inputs=None, variables=None):
        outputsF = [ self.formatInContext(p) for p in outputs ]
        if inputs is not None:
            inputsF  = [ self.formatInContext(p) for p in inputs ]
        else:
            inputsF = None
        if variables is not None:
            variablesF = { k : self.formatInContext(v)
                           for (k, v) in variables.items() }
        else:
            variablesF = None
        buildStmt = BuildStatement(outputsF, rule, inputsF, variablesF)
        self._buildStmts.append(buildStmt)
        return buildStmt

    def write(self, outputNinjaFname="build.ninja"):
        f = open(outputNinjaFname, "w")
        with closing(ninja.Writer(f)) as w:
            w.comment("Variables")
            w.variable("ncpus", "8")
            w.variable("scratchDir", "/scratch")
            if self._grid:
                w.variable("sgeQueue", self._sgeQueue)
                w.variable("grid", "qsub -q $sgeQueue -sync y -cwd -V -b y -e log -o log")
                w.variable("gridSMP", "$grid -pe smp $ncpus")
            else:
                w.variable("grid", "")
                w.variable("gridSMP", "")
            w.newline()
            w.comment("Rules")
            for rule in self._rules.items():
                w.rule(*rule)
                w.newline()
            w.newline()
            w.comment("Build targets")
            for buildStmt in self._buildStmts:
                w.build(buildStmt.outputs, buildStmt.rule, buildStmt.inputs,
                        variables=buildStmt.variables)
                w.newline()

        # Install the bundled resources
        for br in self._resourcesToBundle:
            resSrcPath = getResourcePath(br.resourceName)
            resDestPath = br.destPath or br.resourceName
            mkdirp(op.dirname(resDestPath))
            if not br.substitutions:
                shutil.copy(resSrcPath, resDestPath)
            else:
                tContent = readFile(resSrcPath)
                t = jinja2.Template(tContent)
                writeFile(resDestPath, t.render(br.substitutions))
