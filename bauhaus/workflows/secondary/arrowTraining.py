__all__ = [ "ArrowTrainingWorkflow" ]

from bauhaus import Workflow
from bauhaus.experiment import ResequencingConditionTable

from .mapping import MappingWorkflow


class ArrowTrainingWorkflow(Workflow):
    """
    Run unitem's arrow training, generating a json parameters file
    """
    @staticmethod
    def name():
        return "ArrowTraining"

    @staticmethod
    def conditionTableType():
        return ResequencingConditionTable

    def generate(self, pflow, ct):
        pflow.genRuleOnce("symlink", "ln -s $$(readlink -f $in) $out")
        mapping = MappingWorkflow().generate(pflow, ct)
        ALIGN_XML = "training/{condition}.all_movies.alignmentset.xml"
        REF_FASTA = ALIGN_XML + ".ref.fa"
        REF_FAIDX = REF_FASTA + ".fai"
        inputs = []
        for (condition, alignmentSets) in mapping.items():
            alignmentSet = alignmentSets[0]
            reference = ct.reference(condition)
            with pflow.context("condition", condition):
                inputs.extend(pflow.genBuildStatement([ALIGN_XML], "symlink", [alignmentSet]).outputs)
                inputs.extend(pflow.genBuildStatement([REF_FASTA], "symlink", [reference]).outputs)
                inputs.extend(pflow.genBuildStatement([REF_FAIDX], "symlink", [reference + ".fai"]).outputs)
        pflow.bundleResource("scripts/R/trainArrow.R")
        pflow.genRuleOnce(
            "trainArrow",
            "$grid Rscript --vanilla scripts/R/trainArrow.R training")
        bs = pflow.genBuildStatement(
            [ "fit.json", "Emissions.pdf", "Transitions.pdf" ],
            "trainArrow",
            inputs)
        # return the json parameters file
        return [bs.outputs[0], "trained_condition"]
