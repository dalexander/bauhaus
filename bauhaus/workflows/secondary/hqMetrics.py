from bauhaus import Workflow
from .datasetOps import *
from .mapping import genChunkedMapping
from bauhaus.experiment import (InputType, ResequencingConditionTable)

__all__ = ["HQMetricsWorkflow"]

class HQMetricsWorkflow(Workflow):
    @staticmethod
    def name():
        return "HQMetrics"

    @staticmethod
    def conditionTableType():
        return ResequencingConditionTable

    def generate(self, pflow, ct):
        outputDict = {}
        pflow.bundleScript("R/hqrf_pbbamr.R")
        pflow.genRuleOnce("hqrf_pbbamr",
                          "$grid Rscript --verbose R/hqrf_pbbamr.R $origbam xscraps $in $out")

        pflow.genRuleOnce("bam2bam_nohq",
                          "$grid bam2bam -j $ncpus -b $ncpus --fullHQ --adapters $adapters -o $outPrefix $in")

        pflow.genRuleOnce("subreadset_create",
                          "$grid dataset create $out $in")

        for condition in ct.conditions:
            with pflow.context("condition", condition):
                reference = ct.reference(condition)
                if ct.inputType == InputType.SubreadSet:
                    inputs = ct.inputs(condition)

                    if len(inputs) != 1 :
                        raise NotImplementedError, "condition must be 1:1 with inputs"

                    with pflow.context("movieName", movieName(inputs[0])):
                         # 1. remove HQ regions

                        nohq = pflow.genBuildStatement(
                                [ "{condition}/nohq/{movieName}.subreads.bam" ] ,        # outputs
                                "bam2bam_nohq",  # rule name
                                inputs,
                                dict(outPrefix="{condition}/nohq/{movieName}",
                                     adapters=adaptersFasta(inputs[0]))
                        )

                        nohqsets = pflow.genBuildStatement(
                             [ "{condition}/nohq/{movieName}.subreadset.xml"],
                             "subreadset_create",  # rule name
                             [ "{condition}/nohq/{movieName}.subreads.bam" ] )

                        #print inputs
                        #print nohq.outputs
                        #print nohqsets.outputs
                        # 2. align
#                        noHQAlnSets = genChunkedMapping(pflow, inputs, reference)
                        noHQAlnSets = genChunkedMapping(pflow, nohqsets.outputs, reference)

                        # 3. R script to join the mapping results to the original HQ region start/end markers
                        outputs = [ "{condition}/{movieName}.hqrm.pdf",         # graphical summary output from R script
                                    "{condition}/{movieName}.hqrm_metrics.csv"] # tabular summary output from R script

                        # BEGIN(HACK)
                        # loadRegionsTable() pbbamr method does not support loading a subreadset.xml yet.
                        # but in this case, the source is just a single bam file that has the same root name as
                        # the subreadset.xml. So we'll just hack the file extension for now.
                        origbam = inputs[0] + "" # HACK WARNING FIXME
                        origbam = origbam.replace(".subreadset.xml",".subreads.bam") # HACK WARNING FIXME
                        # END(HACK)

                        b = pflow.genBuildStatement(
                            outputs,        # outputs
                            "hqrf_pbbamr",  # rule name
                            noHQAlnSets,    # inputs
                            dict(origbam=origbam)) # options
                        outputDict[condition] = b.outputs
                else:
                    raise NotImplementedError, "Support not yet implemented for this input type"
        return outputDict
