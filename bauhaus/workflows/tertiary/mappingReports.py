__all__ = [ "MappingReportsWorkflow"]

import os.path as op

from bauhaus import Workflow
from bauhaus.experiment import (InputType, ResequencingConditionTable)
from bauhaus.workflows.secondary import MappingWorkflow
from bauhaus.utils import mkdirp, listConcat

PLOTTING_GROUPS = [ "PbiPlots", "PbiSampledPlots" ]
PLOTTING_TOOLS_ROOT = "/pbi/dept/itg/internaltools" # TODO: ATM we are not bundling these, rather we run them from NFS

def generateConditionsJSON(pflow, ct, alignmentSets):
    """
    Generate a pbcommandR-compatible "condition table" JSON
    """
    pflow.bundleResource("conditions.json",
                         substitutions=dict(
                             conditions=ct.conditions,
                             alignmentSets=alignmentSets,
                             referenceSets={ c : ct.reference(c) for c in ct.conditions }))



def generateToolContract(pflow, toolName):
    """
    Generate a tool contract for running one of the mapping report
    plotting tools from "internaltools"
    """
    destDir = "reports/%s" % toolName
    mkdirp(destDir)
    tcName = "rtc-%s.json" % toolName
    destTcPath = op.join(destDir, tcName)
    noncePath = op.join(destDir, "file.xml")
    conditionsJsonPath = "conditions.json"
    pflow.bundleResource("tool-contracts/%s" % tcName, destTcPath,
                         substitutions = dict(conditions_json=op.abspath(conditionsJsonPath),
                                              output_nonce=op.abspath(noncePath)))
    return noncePath


def generateResequencingPlotToolContracts(pflow):
    for pg in PLOTTING_GROUPS:
        generateToolContract(pflow, pg)

class MappingReportsWorkflow(Workflow):
    """
    Mapping followed by plots/reports
    """
    @staticmethod
    def name():
        return "MappingReports"

    @staticmethod
    def conditionTableType():
        return ResequencingConditionTable

    def generate(self, pflow, ct):
        mapping = MappingWorkflow().generate(pflow, ct)
        alignmentSets = { c : op.abspath(mapping[c][0]) for c in mapping }
        flatMappingOuts = listConcat(list(mapping.values()))
        generateConditionsJSON(pflow, ct, alignmentSets)

        # run plots
        for plotGroup in PLOTTING_GROUPS:
            toolName = plotGroup
            outputNonce = generateToolContract(pflow, toolName)
            pbiPlotsRule = pflow.genRuleOnce(
                toolName,
                "$grid " + PLOTTING_TOOLS_ROOT + "/r/{toolName}.R run-rtc reports/{toolName}/rtc-{toolName}.json".format(toolName=toolName))
            bs = pflow.genBuildStatement(
                [ outputNonce ],
                toolName,
                flatMappingOuts)

        # no output for now
