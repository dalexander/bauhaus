__all__ = [ "InputType",
            "ConditionTable",
            "ResequencingConditionTable",
            "CoverageTitrationConditionTable",
            "UnrolledMappingConditionTable",
            "TableValidationError",
            "InputResolutionError" ]

import csv, os.path as op, sys
import eztable
import io,codecs

# TODO: we are using the non-public method ._get_columns on
# eztable.Table objects.  That method should really be public--maybe
# ask the maintainer?

class TableValidationError(Exception): pass
class InputResolutionError(Exception): pass

class InputType(object):
    SubreadSet            = 1
    AlignmentSet          = 2
    ConsensusReadSet      = 3
    ConsensusAlignmentSet = 4

def _pVariables(tbl):
    return [ k for k in tbl.column_names if k.startswith("p_") ]

def _unique(lst):
    return sorted(set(lst))

class ConditionTable(object):
    """
    Base class for representing Milhouse-style condition tables

    Columns required: Condition, some encoding of input

    There is *no* column encoding secondary protocol; that's the way
    things work in Milhouse but not here.  The "protocol" is specified
    separately from the condition table, which is meant only to encode
    variables and their associated inputs.
    """
    def __init__(self, inputCsv, resolver):
        if not op.isfile(inputCsv):
            raise ValueError("Missing input file: %s" % inputCsv)
        nbytes = min(32, op.getsize(inputCsv))
        raw= open(inputCsv, 'rb').read(nbytes)
        if raw.startswith(codecs.BOM_UTF8):
            raise TableValidationError("Input CSV file is in UTF-8 format. Please convert to ASCII or remove Byte Order Mark (BOM)")
        try:
            with open(inputCsv) as f:
                cr = csv.reader(f)
                allRows = list(cr)
                columnNames, rows = \
                    allRows[0], allRows[1:]
                self.tbl = eztable.Table(columnNames, rows)
        except:
            raise TableValidationError("Input CSV file can't be read/parsed:" + str(sys.exc_info()[0]))
        self._validateTable()
        self._resolveInputs(resolver)

    def _validateTable(self):
        """
        Validate the input CSV file
        Exception on invalid input.

        The intention here is that this can quickly run (by a local
        "validate" rule) to spot errors in the input CSV, before we
        attempt to run the workflow.

        Possible errors that are checked for:
          - Malformed CSV input (can't be parsed)
          - Incorrect header row
          - Too few/many p_ conditions
          - Varying covariate levels under a single condition name
        """
        self._validateConditionsAreHomogeneous()
        self._validateSingleInputEncoding()

    def _validateConditionsAreHomogeneous(self):
        for c in self.conditions:
            condition = self.condition(c)
            for variable in self.variables:
                if len(set(condition._get_column(variable))) != 1:
                    raise TableValidationError(
                        "Conditions must be homogeneous---no variation in variables within a condition")

    def _validateSingleInputEncoding(self):
        cols = self.tbl.column_names
        inputEncodings = 0
        if {"ReportsPath"}.issubset(cols):
            inputEncodings += 1
        if {"RunCode", "ReportsFolder"}.issubset(cols):
            inputEncodings += 1
        if {"SMRTLinkServer", "JobId"}.issubset(cols):
            inputEncodings += 1
        if {"JobPath"}.issubset(cols):
            inputEncodings += 1
        if {"SubreadSet"}.issubset(cols):
            inputEncodings += 1
        if {"AlignmentSet"}.issubset(cols):
            inputEncodings += 1
        if inputEncodings == 0:
          raise TableValidationError("Input data not encoded in condition table. Table requires one and only one column (or pair of columns) as follows: ReportsPath, RunCode+ReportsFolder, SMRTLinkServer+JobId, JobPath, SubreadSet, AlignmentSet")
        if inputEncodings > 1:
            raise TableValidationError("Condition table can only represent input data in one way")

    def _resolveInput(self, resolver, rowRecord):
        cols = self.tbl.column_names
        if {"ReportsPath"}.issubset(cols):
            return resolver.findSubreadSet(rowRecord.ReportsPath)
        elif {"RunCode", "ReportsFolder"}.issubset(cols):
            return resolver.resolveSubreadSet(rowRecord.RunCode, rowRecord.ReportsFolder)
        elif {"SMRTLinkServer", "JobId"}.issubset(cols):
            return resolver.resolveAlignmentSet(rowRecord.SMRTLinkServer, rowRecord.JobId)
        elif {"JobPath"}.issubset(cols):
            return resolver.findAlignmentSet(rowRecord.JobPath)
        elif {"SubreadSet"}.issubset(cols):
            return resolver.ensureSubreadSet(rowRecord.SubreadSet)
        elif {"AlignmentSet"}.issubset(cols):
            return resolver.ensureAlignmentSet(rowRecord.AlignmentSet)


    def _resolveInputs(self, resolver):
        self._inputsByCondition = {}
        for condition in self.conditions:
            subDf = self.condition(condition)
            inputs = []
            for row in subDf:
                inputs.append(self._resolveInput(resolver, row))
            self._inputsByCondition[condition] = inputs

    @property
    def conditions(self):
        return _unique(self.tbl.Condition)

    def condition(self, condition):
        # Get subtable for condition
        return self.tbl.restrict(("Condition",), lambda x: x == condition)

    @property
    def pVariables(self):
        """
        "p variables" encoded in the condition table using column names "p_*"
        """
        return _pVariables(self.tbl)

    @property
    def variables(self):
        """
        "Covariates" are the "p_" variables           -
        """
        return self.pVariables

    def variable(self, condition, variableName):
        """
        Get the value of a variable within a condition

        # TODO: we are using a non-public API method in eztable!
        """
        vals = _unique(self.condition(condition)._get_column(variableName))
        assert len(vals) == 1
        return vals[0]

    @property
    def inputType(self):
        cols = self.tbl.column_names
        if {"ReportsPath"}.issubset(cols) or \
           {"RunCode", "ReportsFolder"}.issubset(cols) or \
           {"SubreadSet"}.issubset(cols):
            return InputType.SubreadSet
        if {"SMRTLinkServer", "JobId"}.issubset(cols) or \
           {"JobPath"}.issubset(cols) or \
           {"AlignmentSet"}.issubset(cols):
            return InputType.AlignmentSet
        raise NotImplementedError, "Input type not recognized/supported"

    def inputs(self, condition):
        return self._inputsByCondition[condition]

    def modelPath(self, condition):
        """
        Use the condition table to locate model parameters files or folders
        """
        if not hasattr(self.condition(condition), 'ModelPath'):
            return None
        paths = _unique(self.condition(condition).ModelPath)
        assert len(paths) == 1
        return paths[0]

    def modelSpec(self, condition):
        """
        Use the condition table to locate model override specification
        """
        if not hasattr(self.condition(condition), "ModelSpec"):
            return None
        specs = _unique(self.condition(condition).ModelSpec)
        assert len(specs) == 1
        return specs[0]

class ResequencingConditionTable(ConditionTable):
    """
    Base class for representing Milhouse-style condition tables for
    resequencing-bases analyses (require a reference, use mapping)
    """
    def _validateGenomeColumnPresent(self):
        if "Genome" not in self.tbl.column_names:
            raise TableValidationError("'Genome' column must be present")

    def _validateTable(self):
        """
        Additional validation: "Genome" column required
        """
        super(ResequencingConditionTable, self)._validateTable()
        self._validateGenomeColumnPresent()

    def _resolveInputs(self, resolver):
        super(ResequencingConditionTable, self)._resolveInputs(resolver)
        self._referenceByCondition = {}
        for condition in self.conditions:
            genome = self.genome(condition)
            self._referenceByCondition[condition] = resolver.resolveReference(genome)

    @property
    def variables(self):
        """
        In addition to the "p_" variables, "Genome" is considered an implicit
        variable in resequencing experiments
        """
        return [ "Genome" ] + self.pVariables

    def genome(self, condition):
        """
        Use the condition table to look up the correct "Genome" based on
        the condition name
        """
        genomes = _unique(self.condition(condition).Genome)
        assert len(genomes) == 1
        return genomes[0]

    def reference(self, condition):
        return self._referenceByCondition[condition]

class CoverageTitrationConditionTable(ResequencingConditionTable):

    def _validateAtLeastOnePVariable(self):
        if len(_pVariables(self.tbl)) < 1:
            raise TableValidationError(
                'There must be at least one covariate ("p_" variable) in the condition table')

    def _validateTable(self):
        super(CoverageTitrationConditionTable, self)._validateTable()
        #self._validateAtLeastOnePVariable()

    def referenceMask(self, condition):
        return self._referenceMaskByCondition[condition]

    def _resolveInputs(self, resolver):
        super(CoverageTitrationConditionTable, self)._resolveInputs(resolver)
        self._referenceMaskByCondition = {}
        for condition in self.conditions:
            genome = self.genome(condition)
            self._referenceMaskByCondition[condition] = resolver.resolveReferenceMask(genome)


class UnrolledMappingConditionTable(ResequencingConditionTable):
    """
    Unrolled mapping requires an unrolled reference.  Ideally we would
    like some kind of flag in the reference information to indicate
    whether this holds or not.  For now, just look for "circular" or
    "unrolled" in the reference name.
    """
    def _validateTable(self):
        super(UnrolledMappingConditionTable, self)._validateTable()
        for genome in self.tbl.Genome:
            if "unrolled" in genome or "circular" in genome:
                continue
            else:
                raise TableValidationError, "Unrolled mapping requires an unrolled reference"
