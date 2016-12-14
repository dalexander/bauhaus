__all__ = [ "VariantCallingWorkflow" ]

from bauhaus import Workflow
from bauhaus.experiment import ResequencingConditionTable

from .mapping import ChunkedMappingWorkflow


# -- Basic variant calling...

def genVariantCalling(pflow, alignedSubreadsSet, modelPath, modelSpec, reference,
                      referenceMask=None, algorithm="arrow", coverageLimit=None):
    """
    Run the variant calling algorithm with all variants filtering
    disabled, producing files:

     - variants{-coverageLimit}.gff
     - consensus{-coverageLimit}.fasta
     - consensus{-coverageLimit}.fastq

    TODO: do we like handling coverageLimit this way, or do we want to leverage PFLow context mechanism?
    """
    if coverageLimit is None:
        CONSENSUS_FASTA       = "{condition}/variant_calling/%s/consensus.fasta"        % algorithm
        CONSENSUS_FASTQ       = "{condition}/variant_calling/%s/consensus.fastq"        % algorithm
        VARIANTS_GFF          = "{condition}/variant_calling/%s/variants.gff"           % algorithm
        MASKED_VARIANTS_GFF   = "{condition}/variant_calling/%s/masked-variants.gff"    % algorithm
        coverageLimitArgument = ""
    else:
        CONSENSUS_FASTA       = "{condition}/variant_calling/%s/consensus-%d.fasta"     % (algorithm, coverageLimit)
        CONSENSUS_FASTQ       = "{condition}/variant_calling/%s/consensus-%d.fastq"     % (algorithm, coverageLimit)
        VARIANTS_GFF          = "{condition}/variant_calling/%s/variants-%d.gff"        % (algorithm, coverageLimit)
        MASKED_VARIANTS_GFF   = "{condition}/variant_calling/%s/masked-variants-%d.gff" % (algorithm, coverageLimit)
        coverageLimitArgument = "-X%d" % coverageLimit

    vcRule = pflow.genRuleOnce(
        "variantCalling",
        "$gridSMP variantCaller $modelPath $modelSpec --algorithm=%s $coverageLimitArgument -x0 -q0 -j $ncpus --reportEffectiveCoverage $in -r $reference -o $out -o $consensusFasta -o $consensusFastq"
        % (algorithm))

    bs = pflow.genBuildStatement(
        [VARIANTS_GFF],
        "variantCalling",
        [alignedSubreadsSet],
        dict(reference=reference,
             modelPath="-P{0}".format(modelPath) if modelPath else "",
             modelSpec="-p{0}".format(modelSpec) if modelSpec else "",
             consensusFasta=CONSENSUS_FASTA,
             consensusFastq=CONSENSUS_FASTQ,
             coverageLimitArgument=coverageLimitArgument))

    if referenceMask is not None:
        maskRule = pflow.genRuleOnce(
            "maskVariantsGff",
            "gffsubtract.pl $in $referenceMask > $out")
        # overwrite bs
        bs = pflow.genBuildStatement(
            [MASKED_VARIANTS_GFF], "maskVariantsGff", [VARIANTS_GFF],
            dict(referenceMask=referenceMask))

    return bs.outputs

def genCoverageSummary(pflow, alignmentSet, reference):
    pflow.genRuleOnce(
        "summarize_coverage",
        "$grid python -m pbreports.report.summarize_coverage.summarize_coverage --region_size=10000 $in $reference $out")
    bs = pflow.genBuildStatement(
        ["{condition}/variant_calling/alignments_summary.gff"],
        "summarize_coverage",
        [alignmentSet],
        dict(reference=reference))
    return bs.outputs


class VariantCallingWorkflow(Workflow):
    """
    Run variant calling on every condition
    """
    @staticmethod
    def name():
        return "VariantCalling"

    @staticmethod
    def conditionTableType():
        return ResequencingConditionTable

    def generate(self, pflow, ct):
        mapping = ChunkedMappingWorkflow().generate(pflow, ct)
        outputDict = {}
        for (condition, alignmentSets) in mapping.iteritems():
            alignmentSet = alignmentSets[0]
            modelPath = ct.modelPath(condition)
            modelSpec = ct.modelSpec(condition)
            reference = ct.reference(condition)
            with pflow.context("condition", condition):
                outputDict[condition] = genVariantCalling(pflow, alignmentSet, modelPath, modelSpec, reference,
                                                          algorithm="arrow")
        return outputDict
