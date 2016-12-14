__all__ = [ "BasicMappingWorkflow", "ChunkedMappingWorkflow"]

import os.path as op

from bauhaus import Workflow
from bauhaus.experiment import (InputType, ResequencingConditionTable)
from .datasetOps import *
from .subreads import genSubreads, genSubreadSetSplit


# Conventions (in the absence of types!):
#
#   - Every generator should returns the list of context-resolved
#     outputs from its final build statement
#


# ----------- Mapping --------------------


# Assumption that should be asserted: each input subreads set comes from a separate movie.

def genMapping(pflow, subreadSets, reference):
    """
    Map the subreads set, without chunking.  This is painfully slow on
    Sequel-scale data.
    """
    mapRule = pflow.genRuleOnce(
        "map",
        "$gridSMP pbalign --tmpDir=$scratchDir --nproc $ncpus $in $reference $out")
    alignmentSets = []
    for subreadSet in subreadSets:
        with pflow.context("entityName", entityName(subreadSet)):
            buildVariables = dict(reference=reference)
            alignmentSets.extend(pflow.genBuildStatement(
                ["{condition}/mapping/{entityName}.alignmentset.xml"],
                "map",
                [subreadSet],
                buildVariables).outputs)
    return genDatasetMergeForCondition(pflow, alignmentSets, "mapping", "alignmentset")

def genMappingCCS(pflow, ccsSets, reference):
    """
    Map... CCS reads.  We really need to make the mapping routine
    generic.  (In the past this actually did require different options
    to blasr/pbalign.  It doesn't anymore.)
    """
    mapRule = pflow.genRuleOnce(
        "mapCCS",
        "$gridSMP pbalign --tmpDir=$scratchDir --nproc $ncpus $in $reference $out")
    alignmentSets = []
    for ccsSet in ccsSets:
        with pflow.context("entityName", entityName(ccsSet)):
            buildVariables = dict(reference=reference)
            alignmentSets.extend(pflow.genBuildStatement(
                ["{condition}/ccs_mapping/{entityName}.consensusalignments.xml"],
                "mapCCS",
                [ccsSet],
                buildVariables).outputs)
    return genDatasetMergeForCondition(pflow, alignmentSets, "ccs_mapping", "consensusalignments")

def genChunkedMapping(pflow, subreadSets, reference, splitFactor=8, doMerge=False, extraBlasrArgs="", extraPbalignArgs=""):
    """
    Break the subreads set into chunks, map the chunks, then
    consolidate the mapped chunks
    """
    extraPbalignArgs = extraPbalignArgs + " --algorithmOptions=\\'%s\\'" % (extraBlasrArgs,) if extraBlasrArgs else ""
    pbalignCmd = "pbalign %s --tmpDir=$scratchDir --nproc $ncpus $in $reference $out" % (extraPbalignArgs,)
    mapRule = pflow.genRuleOnce(
        "map",
        "$gridSMP %s" % (pbalignCmd,))
    alignmentSets = []
    for subreadSet in subreadSets:
        with pflow.context("movieName", movieName(subreadSet)):
            alignmentSetChunks = []
            subreadSetChunks = genSubreadSetSplit(pflow, subreadSet, splitFactor)
            for (i, subreadSetChunk) in enumerate(subreadSetChunks):
                with pflow.context("chunkNum", i):
                    buildVariables = dict(reference=reference)
                    buildStmt = pflow.genBuildStatement(
                        ["{condition}/mapping_chunks/{movieName}.chunk{chunkNum}.alignmentset.xml"],
                        "map",
                        [subreadSetChunk],
                        buildVariables)
                    alignmentSetChunks.extend(buildStmt.outputs)
            if doMerge:
                alignmentSets.extend(
                    genDatasetConsolidateForMovie(pflow, alignmentSetChunks, "mapping", "alignmentset"))
            else:
                alignmentSets.extend(alignmentSetChunks)
    return genDatasetMergeForCondition(pflow, alignmentSets, "mapping", "alignmentset")

# ---------- Workflows -------------
# Unlike the generators above, these guys return a dict of Condition -> resolved outputs
# These really ought to be codified using types.  Damn you Python.
# TODO: do we consider these "tertiary"?

class BasicMappingWorkflow(Workflow):
    """
    Basic mapping---not chunked, just dead simple.
    """
    @staticmethod
    def name():
        return "BasicMapping"

    @staticmethod
    def conditionTableType():
        return ResequencingConditionTable

    def generate(self, pflow, ct):
        outputDict = {}
        for condition in ct.conditions:
            with pflow.context("condition", condition):
                reference = ct.reference(condition)
                if ct.inputType == InputType.SubreadSet:
                    subreadSets = genSubreads(pflow, ct.inputs(condition))
                    outputDict[condition] = genMapping(pflow, subreadSets, reference)
                elif ct.inputType == InputType.AlignmentSet:
                    outputDict[condition] = genAlignmentSetMergeForCondition(pflow, ct.inputs(condition))
                else:
                    raise NotImplementedError, "Support not yet implemented for this input type"
        return outputDict


class ChunkedMappingWorkflow(Workflow):
    """
    Chunked mapping
    """
    @staticmethod
    def name():
        return "ChunkedMapping"

    @staticmethod
    def conditionTableType():
        return ResequencingConditionTable

    def generate(self, pflow, ct):
        outputDict = {}
        for condition in ct.conditions:
            with pflow.context("condition", condition):
                reference = ct.reference(condition)
                if ct.inputType == InputType.SubreadSet:
                    subreadSets = genSubreads(pflow, ct.inputs(condition))
                    outputDict[condition] = genChunkedMapping(pflow, subreadSets, reference, splitFactor=8)
                elif ct.inputType == InputType.AlignmentSet:
                    outputDict[condition] = genDatasetMergeForCondition(pflow, ct.inputs(condition), "mapping", "alignmentset")
                else:
                    raise NotImplementedError, "Support not yet implemented for this input type"
        return outputDict



# -------------------- Demo -------------------

# lambda short insert ecoli runs intended for CCS analysis
if __name__ == '__main__':
    from bauhaus.pflow import PFlow
    inputDataByCondition = { "Replicate1" : ["/pbi/collections/315/3150128/r54008_20160308_001811/1_A01/m54008_160308_002050.subreadset.xml"],
                             "Replicate2" : ["/pbi/collections/315/3150128/r54008_20160308_001811/2_B01/m54008_160308_053311.subreadset.xml"] }

    reference = "/mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta"

    pflow = PFlow()
    for (condition, inputSubreadSets) in inputDataByCondition.iteritems():
        with pflow.context("condition", condition):
            #genMapping(pflow, inputSubreadSets, reference)
            genChunkedMapping(pflow, inputSubreadSets, reference, 8)
    pflow.write("build.ninja")
