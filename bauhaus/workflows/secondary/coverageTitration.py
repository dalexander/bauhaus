__all__ = [ "CoverageTitrationWorkflow",
            "CoverageTitrationReportsWorkflow" ]

from bauhaus import Workflow
from bauhaus.utils import listConcat
from bauhaus.experiment import CoverageTitrationConditionTable

from .variantCalling import genVariantCalling, genCoverageSummary
from .mapping import MappingWorkflow

from collections import defaultdict

COVERAGE_LEVELS = [5, 10, 15, 20, 30, 40, 50, 60, 80, 100]


class CoverageTitrationWorkflow(Workflow):
    """
    Run coverage titration on every condition
    """
    @staticmethod
    def name():
        return "CoverageTitration"

    @staticmethod
    def conditionTableType():
        return CoverageTitrationConditionTable

    def generate(self, pflow, ct):
        mapping = MappingWorkflow().generate(pflow, ct)
        outputDict = defaultdict(list)
        for (condition, alignmentSets) in mapping.iteritems():
            alignmentSet = alignmentSets[0]
            modelPath = ct.modelPath(condition)
            modelSpec = ct.modelSpec(condition)
            reference = ct.reference(condition)
            referenceMask = ct.referenceMask(condition)
            with pflow.context("condition", condition):
                coverageSummary = genCoverageSummary(pflow, alignmentSet, reference)[0]
                outputDict[condition].append(coverageSummary)
                algorithm = "arrow" # TODO: make configurable!
                with pflow.context("consensusAlgorithm", algorithm):
                    for x in COVERAGE_LEVELS:
                        outputDict[condition].extend(
                            genVariantCalling(pflow, alignmentSet, modelPath, modelSpec, reference,
                                              referenceMask=referenceMask,
                                              algorithm=algorithm,
                                              coverageLimit=x))
        return outputDict


class CoverageTitrationReportsWorkflow(Workflow):
    """
    Run coverage titration on every condition, then generate
    comparative analyses/plots of the different conditions
    """
    @staticmethod
    def name():
        return "CoverageTitrationReports"

    @staticmethod
    def conditionTableType():
        return CoverageTitrationConditionTable

    def generate(self, pflow, ct):
        pflow.bundleResource("scripts/R/coverageTitrationPlots.R")
        ctOuts = CoverageTitrationWorkflow().generate(pflow, ct)
        flatCtOuts = listConcat(ctOuts.values())
        ctSummaryRule = pflow.genRuleOnce(
            "coverageTitrationSummaryAnalysis",
            "Rscript --vanilla scripts/R/coverageTitrationPlots.R .")
        bs = pflow.genBuildStatement(
            [ "coverage-titration.csv", "coverage-titration.pdf"],
            "coverageTitrationSummaryAnalysis",
            flatCtOuts)
