__all__ = [ "genSubreads", "genSubreadSetSplit" ]

import os.path as op, re

from bauhaus import Workflow
from .datasetOps import *
from bauhaus.experiment import (InputType, ConditionTable)

def genSubreads(pflow, remoteSubreadSets):
    """
    Copy subreads to the workflow dir, while fixing relative paths
    """
    copyRule = pflow.genRuleOnce(
        "copySubreadsDataset",
        "$grid dataset create $out $in")
    localSubreadSets = []
    for rss in remoteSubreadSets:
        bn = op.basename(rss)
        if bn.endswith(".subreads.bam"):
            bn = re.sub(".subreads.bam$", ".subreadset.xml", bn)
        localName = "{condition}/subreads/%s" % bn
        subreadCopyStmt = pflow.genBuildStatement(
            ["{condition}/subreads/%s" % bn],
            "copySubreadsDataset",
            [rss])
        localSubreadSets.extend(subreadCopyStmt.outputs)
    return localSubreadSets


def genSubreadSetSplit(pflow, subreadSet, splitFactor):
    # Split by ZMWs.  Returns: [dset]
    assert splitFactor >= 1
    pflow.genRuleOnce(
            "splitByZmw",
            "$grid dataset split --zmws --targetSize 1 --chunks %d --outdir $outdir $in" % (splitFactor,))
    movie = movieName(subreadSet)
    splitOutputs =  [ "{condition}/subreads_chunks/%s.chunk%d.subreadset.xml" % (movie, i)
                      for i in xrange(splitFactor) ]
    buildStmt = pflow.genBuildStatement(splitOutputs,
                                        "splitByZmw",
                                        [subreadSet],
                                        variables={"outdir": "{condition}/subreads_chunks"})
    return buildStmt.outputs
