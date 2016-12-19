import os.path as op
import re
from xml.etree import ElementTree

def movieName(filename):
    """
    The movieName is reckoned as the first .-delimited part of the
    basename.

    (Filename is:

     dir/movieName[.chunkName].dsettype.xml    )
    """

    base = op.basename(filename)
    fields = base.split(".")
    assert len(fields) > 2 and fields[-1] == "xml"
    if len(fields) > 3 and re.match(r"chunk[0-9]+", fields[-3]):
        movieName = ".".join(fields[:-3])
    else:
        movieName = ".".join(fields[:-2])
    return movieName

def entityName(filename):
    """
    An entity is either a movie or a part (chunk) of a movie

    This is encoded in the filename as:

     dir/movieName[.chunkName].dsettype.xml
    """
    base = op.basename(filename)
    fields = base.split(".")
    assert len(fields) in (3, 4) and fields[-1] == "xml"
    return ".".join(fields[:-2])

def subreadsBam(subreadSet):
    et = ElementTree.parse(subreadSet)
    q = et.findall('.//{http://pacificbiosciences.com/PacBioBaseDataModel.xsd}ExternalResource[@MetaType="PacBio.SubreadFile.SubreadBamFile"]')
    assert len(q) == 1
    subreadsBamFileName=q[0].attrib["ResourceId"]
    if not op.isabs(subreadsBamFileName):
        subreadsBamFileName = op.join(op.dirname(subreadSet), subreadsBamFileName)
    return subreadsBamFileName

def scrapsBam(subreadSet):
    et = ElementTree.parse(subreadSet)
    q = et.findall('.//{http://pacificbiosciences.com/PacBioBaseDataModel.xsd}ExternalResource[@MetaType="PacBio.SubreadFile.ScrapsBamFile"]')
    if len(q) == 0:
        raise IOError,"scraps.bam not found in " + subreadSet
    if len(q) > 1:
        print "While parsing " + subreadSet + " found " + str(len(q)) + " scraps.bam resources:"
        for x in q:
            print x.attrib["ResourceId"]
        raise IOError,"too many scraps.bams found"
    assert len(q) == 1
    scrapsBamFileName=q[0].attrib["ResourceId"]
    if not op.isabs(scrapsBamFileName):
        scrapsBamFileName = op.join(op.dirname(subreadSet), scrapsBamFileName)
    return scrapsBamFileName

def reportsDirectory(subreadSet):
    subreadsBamFileName = subreadsBam(subreadSet)
    assert op.isabs(subreadsBamFileName)
    return op.dirname(subreadsBamFileName)

def adaptersFasta(subreadSet):
    """
    Presently we poke around in the run dir to find the `adapters.fasta`.
    In the future, this will be indicated as an external resource in the dataset XML.
    """
    return op.join(reportsDirectory(subreadSet), movieName(subreadSet) + ".adapters.fasta")


def extractedAdaptersFasta(subreadSet):
    """
    Look in the XML file for the adapters.fasta file. If not found, look in the same
    folder as the XML file.
    """

    et = ElementTree.parse(subreadSet)
    q = et.findall('.//{http://pacificbiosciences.com/PacBioBaseDataModel.xsd}ExternalResource[@MetaType="PacBio.SubreadFile.AdapterFastaFile"]')
    if len(q) == 0:
        fastaFile = op.join(reportsDirectory(subreadSet), movieName(subreadSet) + ".adapters.fasta")
    elif len(q) > 1:
        print "While parsing " + subreadSet + " found " + str(len(q)) + " adapters resources:"
        for x in q:
            print x.attrib["ResourceId"]
        raise IOError,"too many adapsters.fasta found"
    else:
        fastaFile=q[0].attrib["ResourceId"]
    if not op.isabs(fastaFile):
        fastaFile = op.join(op.dirname(subreadSet), fastaFile)
    return fastaFile

# ------- Split/merge -------------


# Conventions (in the absence of types!):
#
#   - Every generator should returns the list of context-resolved
#     outputs from its final build statement
#


# ---------

def genDatasetMergeForCondition(pflow, inputDatasets, taskTypeName, datasetTypeName):
    # Just merge datasets, for condition
    pflow.genRuleOnce("mergeDatasetsForCondition",
                      "$grid dataset merge $out $in")
    outputs = [ "{condition}/%s/all_movies.%s.xml" % (taskTypeName, datasetTypeName) ]
    buildStmtMerge = pflow.genBuildStatement(
        outputs,
        "mergeDatasetsForCondition",
        inputDatasets)
    return buildStmtMerge.outputs

def genDatasetConsolidateForMovie(pflow, datasets, taskTypeName, datasetTypeName):
    # 1) Merge datasets and 2) physically consolidate BAM files
    # 1.
    pflow.genRuleOnce("mergeDatasetsForMovie",
                      "$grid dataset merge $out $in")
    outputs = [ "{condition}/%s/{movieName}_preconsolidate.%s.xml" % \
                (taskTypeName, datasetTypeName) ]
    buildStmtMerge = pflow.genBuildStatement(
        outputs,
        "mergeDatasetsForMovie",
        datasets)
    mergedDataset = buildStmtMerge.outputs[0]
    # 2.
    pflow.genRuleOnce(
        "consolidateDatasetsForMovie",
        "$grid dataset consolidate $in $outBam $out")
    consolidatedDataset = "{condition}/%s/{movieName}.%s.xml" % (taskTypeName, datasetTypeName)
    # FIXME: this is gross.  Let's do something better.  For example,
    # we could have dataset subclasses that know how to perform
    # operations.
    if datasetTypeName == "alignmentset":
        bamTypeName = "aligned_subreads"
    elif datasetTypeName == "unrolledalignmentset":
        bamTypeName = "aligned_unrolledreads"
    elif datasetTypeName == "consensusreadset":
        bamTypeName = "ccs"
    else:
        raise Exception, "Unsupported dataset type..."
    consolidatedBam = "{condition}/%s/{movieName}.%s.bam" % (taskTypeName, bamTypeName)
    consolidatedDataset = "{condition}/%s/{movieName}.%s.xml" % (taskTypeName, datasetTypeName)
    buildStmtConsolidate = pflow.genBuildStatement(
        [consolidatedDataset],
        "consolidateDatasetsForMovie",
        [mergedDataset],
        dict(outBam=consolidatedBam))
    return buildStmtConsolidate.outputs
