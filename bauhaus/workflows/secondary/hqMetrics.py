from bauhaus import Workflow
from bauhaus.workflows.secondary.datasetOps import *
from bauhaus.workflows.secondary.mapping import genChunkedMapping
from bauhaus.experiment import (InputType, ResequencingConditionTable)
import os.path

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
        rscript = "R/hqrf_pbbamr.R"
        pflow.bundleScript(rscript)
        controls = "fasta/600bp_Control_c2.fasta"
        pflow.bundleScript(controls)

        # R script merges  nohq_aligned reads with HQR regions
        pflow.genRuleOnce("hqrf_pbbamr",
                          "$grid Rscript --verbose " + rscript + " $subreadsbam $scrapsbam $in $out")

        # strips out HQR region information (forming zmw reads) and re-does adapter finding on the zmw reads
        pflow.genRuleOnce("bam2bam_nohq",
                          "$grid bam2bam -j $ncpus -b $ncpus --fullHQ --adapters $adapters -o $outPrefix $in")

        # converts BAZ file to subreads.bam, doing adapter finding and control filtering
        pflow.genRuleOnce("baz2bam",
                          "$grid baz2bam -j $ncpus -b $ncpus --minSnr=0 --metadata $metadata --adapters $adapters --controls $controls -o $outPrefix $in")

        # runs basecaller on trace file, outputting BAZ file
        pflow.genRuleOnce("t2b",
                          "$grid stdbuf -o 0 $basecaller --internal --input $in --output $out --logoutput $log")

        # creates a subreadset dataset from bam files
        pflow.genRuleOnce("subreadset_create",
                          "$grid dataset create $out $in")

        for condition in ct.conditions:
            with pflow.context("condition", condition):
                reference = ct.reference(condition)

                inputs = ct.inputs(condition)

                # the workflow can start with TraceH5File or SubreadSet
                # TODO : add BAZFile as a starting point
                if ct.inputType == InputType.TraceH5File:
                    if len(inputs) != 1 :
                        raise NotImplementedError, "condition must be 1:1 with inputs"

                    trch5 = inputs[0] + ""
                    (d,b) = os.path.split(trch5)
                    movieName_ = b.replace(".trc.h5","")

                    # TODO: find meta data and/or adapters in a better way (not assume a naming convention)
                    metadata = os.path.join(d,"." + movieName_ + ".run.metadata.xml")
                    adaptersFile = os.path.join(d,movieName_ + ".adapters.fasta")
                    with pflow.context("movieName", movieName_ ):

                        baz = "{condition}/t2b/{movieName}.baz"
                        pflow.genBuildStatement(
                            [ baz ] ,        # outputs
                            "t2b",
                            [ trch5 ],
                            dict(basecaller="basecaller-console-app",
                                 log="{condition}/t2b/basecaller-console-app.log")
                        )

                        subreadset = "{condition}/hq/{movieName}.subreadset.xml"
                        pflow.genBuildStatement(
                            [ subreadset  ] ,        # outputs
                            "baz2bam",
                            [ baz ],
                            dict(metadata=metadata,
                                 controls=controls,
                                 adapters=adaptersFile,
                                 outPrefix="{condition}/hq/{movieName}")
                        )

                elif ct.inputType == InputType.SubreadSet:
                    inputs = ct.inputs(condition)

                    if len(inputs) != 1 :
                        raise NotImplementedError, "condition must be 1:1 with inputs"

                    subreadset = inputs[0] + ""
                    movieName_ = movieName(subreadset)
                    adaptersFile = extractedAdaptersFasta(subreadset)

                else:
                    raise NotImplementedError, "Support not yet implemented for this input type %s" % inputs[0]


                with pflow.context("movieName", movieName_ ):

                    # 1. remove HQ regions
                    nohq = pflow.genBuildStatement(
                            [ "{condition}/nohq/{movieName}.subreads.bam" ] ,        # outputs
                            "bam2bam_nohq",  # rule name
                            [ subreadset ],
                            dict(outPrefix="{condition}/nohq/{movieName}",
                                 adapters=adaptersFile)
                    )

                    nohqsets = pflow.genBuildStatement(
                         [ "{condition}/nohq/{movieName}.subreadset.xml"],
                         "subreadset_create",  # rule name
                         [ "{condition}/nohq/{movieName}.subreads.bam" ] )

                    # 2. align
                    noHQAlnSets = genChunkedMapping(pflow, nohqsets.outputs, reference)

                    # 3. R script to join the mapping results to the original HQ region start/end markers
                    outputs = [ "{condition}/{movieName}.hqrm.pdf",         # graphical summary output from R script
                                "{condition}/{movieName}.hqrm_metrics.csv"] # tabular summary output from R script

                    # BEGIN(HACK)
                    # This is needed because pbbamr doesn't support subreadset.xml yet.
                    # If the subreadset is in the condition table, then it can be parsed at bauhaus time (i.e. before run time).
                    # but if the subreadset is generated dynamically, then it can't be parsed at bauhaus time,
                    # so we guess as to what the subreads.bam and scraps.bam file paths will be.
                    if ct.inputType == InputType.SubreadSet:
                        subreadsBamFile = subreadsBam(subreadset)
                        scrapsBamFile   = scrapsBam(subreadset)
                    else:
                        subreadsBamFile = subreadset.replace(".subreadset.xml",".subreads.bam") # HACK WARNING FIXME
                        scrapsBamFile   = subreadset.replace(".subreadset.xml",".scraps.bam") # HACK WARNING FIXME
                    # END(HACK)

                    b = pflow.genBuildStatement(
                        outputs,        # outputs
                        "hqrf_pbbamr",  # rule name
                        noHQAlnSets,    # inputs
                        dict(subreadsbam=subreadsBamFile,  # options
                             scrapsbam=scrapsBamFile))

                    outputDict[condition] = b.outputs
        return outputDict


if __name__ == "__main__":
    wf    = HQMetricsWorkflow()
    pflow = PFlow()
    ct    = ResequencingConditionTable()
    wf.generate(pflow, ct)
