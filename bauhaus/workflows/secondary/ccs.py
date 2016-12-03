__all__ = [ "BasicCCSWorkflow",
            "ChunkedCCSWorkflow",
            "CCSMappingWorkflow",
            "CCSMappingReportsWorkflow" ]

import os.path as op

from bauhaus.experiment import InputType, ConditionTable, ResequencingConditionTable
from bauhaus import Workflow
from bauhaus.utils import listConcat

from .datasetOps import *
from .mapping import genMappingCCS
from .subreads import genSubreads, genSubreadSetSplit

def formatModelParams(buildVariables, modelPath, modelSpec):
    buildVariables["modelPath"] = "--modelPath={0}".format(modelPath) if modelPath else ""
    buildVariables["modelSpec"] = "--modelSpec={0}".format(modelSpec) if modelSpec else ""
    return buildVariables

def genCCS(pflow, subreadSets, modelPath, modelSpec):
    ccsRule = pflow.genRuleOnce(
        "ccs",
        "$gridSMP ccs --force --numThreads=$ncpus --reportFile=$ccsDiagnostics $modelPath $modelSpec $in $outBam && " +
        "dataset create --type ConsensusReadSet $out $outBam")
    ccsSets = []
    for subreadSet in subreadSets:
        with pflow.context("movieName", movieName(subreadSet)):
            outBam = pflow.formatInContext("{condition}/ccs/{movieName}.ccs.bam")
            ccsDiagnostics= pflow.formatInContext("{condition}/{movieName}.ccs-report.txt")
            buildVariables = dict(outBam=outBam, ccsDiagnostics=ccsDiagnostics)
            buildVariables = formatModelParams(buildVariables, modelPath, modelSpec)
            ccsSets.extend(pflow.genBuildStatement(
                ["{condition}/ccs/{movieName}.consensusreadset.xml"],
                "ccs",
                [subreadSet],
                buildVariables))
    return ccsSets

def genChunkedCCS(pflow, subreadSets, modelPath, modelSpec, splitFactor=8, doMerge=True):
    ccsRule = pflow.genRuleOnce(
        "ccs",
        "$gridSMP ccs --force --numThreads=$ncpus --reportFile=$ccsDiagnostics $modelPath $modelSpec $in $outBam && " +
        "dataset create --type ConsensusReadSet $out $outBam")
    ccsSets = []
    for subreadSet in subreadSets:
        with pflow.context("movieName", movieName(subreadSet)):
            ccsSetChunks = []
            subreadSetChunks = genSubreadSetSplit(pflow, subreadSet, splitFactor)
            for (i, subreadSetChunk) in enumerate(subreadSetChunks):
                with pflow.context("chunkNum", i):
                    outBam = pflow.formatInContext("{condition}/ccs_chunks/{movieName}.chunk{chunkNum}.ccs.bam")
                    ccsDiagnostics = pflow.formatInContext("{condition}/ccs_chunks/{movieName}.chunk{chunkNum}.ccs-report.txt")
                    buildVariables = dict(outBam=outBam, ccsDiagnostics=ccsDiagnostics)
                    buildVariables = formatModelParams(buildVariables, modelPath, modelSpec)
                    buildStmt = pflow.genBuildStatement(
                        ["{condition}/ccs_chunks/{movieName}.chunk{chunkNum}.consensusreadset.xml"],
                        "ccs",
                        [subreadSetChunk],
                        buildVariables)
                    ccsSetChunks.extend(buildStmt.outputs)
            # Consolidate or merge the CCS chunks from this movie
            if doMerge:
                ccsSets.extend(
                    genDatasetConsolidateForMovie(pflow, ccsSetChunks, "ccs", "consensusreadset"))
            else:
                ccsSets.extend(ccsSetChunks)
    # Merge for the condition
    if doMerge:
        return genDatasetMergeForCondition(pflow, ccsSets, "ccs", "consensusreadset")
    else:
        return ccsSets


def genCCSAndMapping(pflow, subreadSets, modelPath, modelSpec, reference):
    ccsSets = genChunkedCCS(pflow, subreadSets, modelPath, modelSpec, doMerge=False)
    alignmentSets = genMappingCCS(pflow, ccsSets, reference)
    return alignmentSets


# -- CCS
class BasicCCSWorkflow(Workflow):
    """
    Basic CCS---not chunked, just dead simple.
    """
    @staticmethod
    def name():
        return "BasicCCS"

    @staticmethod
    def conditionTableType():
        return ConditionTable

    def generate(self, pflow, ct):
        outputDict = {}
        for condition in ct.conditions:
            with pflow.context("condition", condition):
                assert ct.inputType == InputType.SubreadSet # TODO: move to validation.
                subreadSets = genSubreads(pflow, ct.inputs(condition))
                modelPath = ct.modelPath(condition)
                modelSpec = ct.modelSpec(condition)
                outputDict[condition] = genCCS(pflow, subreadSets, modelPath, modelSpec)
        return outputDict


class ChunkedCCSWorkflow(Workflow):
    """
    Chunked CCS---split up into multiple jobs
    """
    @staticmethod
    def name():
        return "ChunkedCCS"

    @staticmethod
    def conditionTableType():
        return ConditionTable

    def generate(self, pflow, ct):
        outputDict = {}
        for condition in ct.conditions:
            with pflow.context("condition", condition):
                assert ct.inputType == InputType.SubreadSet # TODO: move to validation.
                subreadSets = genSubreads(pflow, ct.inputs(condition))
                modelPath = ct.modelPath(condition)
                modelSpec = ct.modelSpec(condition)
                outputDict[condition] = genChunkedCCS(pflow, subreadSets, modelPath, modelSpec)
        return outputDict


# -- CCS mapping
class CCSMappingWorkflow(Workflow):
    """
    CCS + mapping
    """
    @staticmethod
    def name():
        return "CCSMapping"

    @staticmethod
    def conditionTableType():
        return ResequencingConditionTable

    def generate(self, pflow, ct):
        outputDict = {}
        for condition in ct.conditions:
            with pflow.context("condition", condition):
                assert ct.inputType == InputType.SubreadSet # TODO: move to validation.
                subreadSets = genSubreads(pflow, ct.inputs(condition))
                modelPath = ct.modelPath(condition)
                modelSpec = ct.modelSpec(condition)
                reference = ct.reference(condition)
                outputDict[condition] = genCCSAndMapping(pflow, subreadSets, modelPath, modelSpec, reference)
        return outputDict



# -- CCS mapping + reports
class CCSMappingReportsWorkflow(Workflow):

    @staticmethod
    def name():
        return "CCSMappingReports"

    @staticmethod
    def conditionTableType():
        return ResequencingConditionTable

    def generate(self, pflow, ct):
        pflow.bundleScript("R/ccsMappingPlots.R")
        ccsMappingOutputs = CCSMappingWorkflow().generate(pflow, ct)
        flatOutputs = listConcat(ccsMappingOutputs.values())

        ccsSummaryRule = pflow.genRuleOnce(
            "ccsMappingSummaryAnalysis",
            "Rscript --vanilla R/ccsMappingPlots.R .")
        bs = pflow.genBuildStatement(
            [ "ccs-mapping.pdf"],
            "ccsMappingSummaryAnalysis",
            flatOutputs)
