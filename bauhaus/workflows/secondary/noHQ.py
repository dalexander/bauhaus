from bauhaus import Workflow
from .datasetOps import *
from .mapping import genChunkedMapping
from bauhaus.experiment import (InputType, ResequencingConditionTable)

__all__ = ["NoHQMappingWorkflow"]


def genNoHQSubreads(pflow, subreadSets):
    """
    Use bam2bam to recompute the subreads, as if the HQRF had been
    disabled.
    """
    bam2bamNoHQRule = pflow.genRuleOnce(
        "deHQify",
        "$gridSMP bam2bam -b $ncpus -j $ncpus --fullHQ --adapters $adapters $subreadsIn $scrapsIn -o $outPrefix &&" +
        "dataset create --type SubreadSet $out $outBam")
    noHQSubreadSets = []
    for subreadSet in subreadSets:
        with pflow.context("movieName", movieName(subreadSet)):
            adapters= adaptersFasta(subreadSet)
            noHQSubreadSets.extend(
                pflow.genBuildStatement(
                    ["{condition}/subreads_noHQ/{movieName}.subreadset.xml"],
                    "deHQify",
                    [subreadSet],
                    dict(outBam="{condition}/subreads_noHQ/{movieName}.subreads.bam",
                         outPrefix="{condition}/subreads_noHQ/{movieName}",
                         subreadsIn=subreadsBam(subreadSet),
                         scrapsIn=scrapsBam(subreadSet),
                         adapters=adapters)
                ).outputs)
    return noHQSubreadSets


class NoHQMappingWorkflow(Workflow):
    @staticmethod
    def name():
        return "NoHQMapping"

    @staticmethod
    def conditionTableType():
        return ResequencingConditionTable

    def generate(self, pflow, ct):
        outputDict = {}
        for condition in ct.conditions:
            with pflow.context("condition", condition):
                reference = ct.reference(condition)
                if ct.inputType == InputType.SubreadSet:
                    noHQSubreads = genNoHQSubreads(pflow, ct.inputs(condition))
                    noHQAlnSets = genChunkedMapping(pflow, noHQSubreads, reference)
                    outputDict[condition] = noHQAlnSets
                else:
                    raise NotImplementedError("Support not yet implemented for this input type")
        return outputDict
