from bauhaus import Workflow
from bauhaus.experiment import (InputType, UnrolledMappingConditionTable)
from .mapping import genChunkedMapping

from .datasetOps import *

__all__ = ["UnrolledNoHQMappingWorkflow"]

# --

def genUnrolledNoHQMapping(pflow, subreadSets, reference, splitFactor=8, doMerge=False):
    unrolledPbalignOptions = "--noSplitSubreads --hitPolicy=leftmost"
    unrolledBlasrOptions = "--bestn 1 --forwardOnly --fastMaxInterval --maxAnchorsPerPosition 30000 --ignoreHQRegions --minPctIdentity 60"
    return genChunkedMapping(pflow, subreadSets, reference, splitFactor, doMerge,
                             extraPbalignArgs=unrolledPbalignOptions,
                             extraBlasrArgs=unrolledBlasrOptions)

# --
# TODO: I am actually not even sure whether BLASR supports unrolled
# mapping while still respecting the HQ region.  We would like to have
# such a workflow...


class UnrolledNoHQMappingWorkflow(Workflow):
    @staticmethod
    def name():
        return "UnrolledNoHQMapping"

    @staticmethod
    def conditionTableType():
        return UnrolledMappingConditionTable

    def generate(self, pflow, ct):
        outputDict = {}
        for condition in ct.conditions:
            with pflow.context("condition", condition):
                reference = ct.reference(condition)
                if ct.inputType == InputType.SubreadSet:
                    outputDict[condition] = genUnrolledNoHQMapping(
                        pflow, ct.inputs(condition), reference, splitFactor=8)
                else:
                    raise NotImplementedError, "Support not yet implemented for this input type"
        return outputDict
