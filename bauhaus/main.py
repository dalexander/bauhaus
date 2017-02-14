from __future__ import print_function
import argparse, shutil, sys, os, os.path as op

from bauhaus.pbls2 import Resolver, MockResolver
from bauhaus.pflow import PFlow
from bauhaus.workflows import availableWorkflows
from bauhaus.utils import mkdirp
from bauhaus.experiment import (InputResolutionError, TableValidationError)
from bauhaus.__version__ import __VERSION__

def doWorkflowHelp(args):
    if not args.workflow :
        print("Please pick a workflow with the --workflow option:")
        print(' '.join(list(availableWorkflows.keys())))
        return

    wfg = availableWorkflows[args.workflow]()
    wfg.help()

def doValidate(args):
    if args.mockResolver:
        r = MockResolver()
    else:
        r = Resolver()
    wfg = availableWorkflows[args.workflow]()
    ct = wfg.conditionTableType()(args.conditionTable, r)
    return wfg, ct

def doGenerate(args, wfg, ct):
    pflow = PFlow()
    if args.noGrid:
        pflow.noGrid()
    else:
        pflow.grid(args.sgeQueue)
    pflow.chunks = args.chunks
    wfg.generate(pflow, ct)
    pflow.write("build.ninja")
    print('Runnable workflow written to directory "%s"' % args.outputDirectory)

def doRun(args):
    raise NotImplementedError

def parseArgs():
    parser = argparse.ArgumentParser(prog="bauhaus")
    parser.add_argument("--version", action="version", version=__VERSION__)
    parser.add_argument(
        "--conditionTable", "-t",
        action="store", metavar="CONDITION_TABLE.CSV",
        required=True,
        type=op.abspath)
    parser.add_argument(
        "--workflow", "-w",
        action="store", type=str,
        required=True,
        choices = list(availableWorkflows.keys()))
    parser.add_argument(
        "--mockResolver", "-m",
        action="store_true",
        help="Use mock pbls2 resolver (for testing purposes)")
    parser.add_argument(
        "--pdb", action="store_true",
        help="Drop into debugger on exception")
    parser.add_argument(
        "--outputDirectory", "-o",
        default="out",
        action="store", type=str)
    parser.add_argument(
        "--noGrid", action="store_true",
        help="Disable the qsub submission to the grid")

    parser.add_argument(
        "-q", "--sgeQueue",
        default="default", type=str,
        help="Specify destination SGE queue for workflow tasks")
    parser.add_argument(
        "--chunks", type=int, default=8,
        help="The number of chunks that should be used for scatter-gather compatible workflows (0 means disable chunking)")

    subparsers = parser.add_subparsers(help="sub-command help", dest="command")
    subparsers.add_parser("help", help="Help for the given work flow")
    subparsers.add_parser("validate", help="Validate the condition table")
    subparsers.add_parser("generate", help="Generate the ninja script to run the workflow")
    subparsers.add_parser("run", help="Run the workflow")

    args = parser.parse_args()
    return args


def _main(args):
    #print args

    if args.command == "help":
        doWorkflowHelp(args)
        return

    if args.command in ("validate", "generate", "run"):
        try:
            wfg, ct = doValidate(args)
            print("Validation and input resolution succeeded.")
            if args.command == "validate": return 0
        except (TableValidationError, InputResolutionError) as e:
            print("Condition table validation error:", e)
            return 1

    # Set up workflow directory

    mkdirp(args.outputDirectory)
    shutil.copyfile(args.conditionTable, op.join(args.outputDirectory, "condition-table.csv"))
    os.chdir(args.outputDirectory)
    mkdirp("log")


    if args.command == "generate":
        doGenerate(args, wfg, ct)
    elif args.command == "run":
        doRun(args)


def main():
    args = parseArgs()
    if args.pdb:
        try:
            import ipdb
            with ipdb.launch_ipdb_on_exception():
                _main(args)
            return 0
        except ImportError:
            return _main(args)
    else:
        return _main(args)


if __name__ == '__main__':
    sys.exit(main())
