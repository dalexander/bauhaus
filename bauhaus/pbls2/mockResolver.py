__all__ = [ "MockResolver" ]

import os.path as op
from .resolver import _isRuncode
from .exceptions import *

class MockResolver(object):
    # For testing purposes

    REFERENCE_MASKS_ROOT = "/mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks"
    REFERENCES_ROOT = "/mnt/secondary/iSmrtanalysis/current/common/references"

    def __init__(self):
        pass

    def resolveSubreadSet(self, runCode, reportsFolder=""):
        if not _isRuncode(runCode):
            raise ValueError('Argument "%s" does not appear to be a runcode' % runCode)
        lookup = \
            { "3150128-0001" : "/pbi/collections/315/3150128/r54008_20160308_001811/1_A01/m54008_160308_002050.subreadset.xml" ,
              "3150128-0002" : "/pbi/collections/315/3150128/r54008_20160308_001811/2_B01/m54008_160308_053311.subreadset.xml" ,
              "3150122-0001" : "/pbi/collections/315/3150122/r54011_20160305_235615/1_A01/m54011_160305_235923.subreadset.xml" ,
              "3150122-0002" : "/pbi/collections/315/3150122/r54011_20160305_235615/2_B01/m54011_160306_050740.subreadset.xml" }
        if runCode not in lookup or reportsFolder != "":
            raise DataNotFound(runCode)
        return lookup[runCode]

    def resolveReference(self, referenceName):
        if referenceName not in ["lambdaNEB", "ecoliK12_pbi_March2013"]:
            raise DataNotFound("Reference not found: %s" % referenceName)
        referenceFasta = op.join(self.REFERENCES_ROOT, referenceName, "sequence", referenceName + ".fasta")
        return referenceFasta

    def resolveReferenceMask(self, referenceName):
        if referenceName not in ["lambdaNEB", "ecoliK12_pbi_March2013"]:
            raise DataNotFound("Reference mask not found: %s" % referenceName)
        return op.join(self.REFERENCE_MASKS_ROOT, referenceName + "-mask.gff")

    def resolveJob(self, smrtLinkServer, jobId):
        lookup = { ("smrtlink-beta", "4110") : "/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite/userdata/jobs_root/004/004110",
                   ("smrtlink-beta", "4111") : "/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite/userdata/jobs_root/004/004111",
                   ("smrtlink-beta", "4183") : "/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite/userdata/jobs_root/004/004183",
                   ("smrtlink-beta", "4206") : "/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite/userdata/jobs_root/004/004206" }
        if (smrtLinkServer, jobId) not in lookup:
            raise DataNotFound("Job not found: %s:%s" % (smrtLinkServer, jobId))
        else:
            return lookup[(smrtLinkServer, jobId)]

    def resolveAlignmentSet(self, smrtLinkServer, jobId):
        jobDir = self.resolveJob(smrtLinkServer, jobId)
        return op.join(jobDir, "tasks/pbalign.tasks.consolidate_bam-0/final.alignmentset.alignmentset.xml")

    def ensureSubreadSet(self, subreadSet):
        if not subreadSet.endswith(".subreadset.xml"):
            raise InvalidDataset, "%s not a subreadset"
        else:
            return subreadSet

    def ensureAlignmentSet(self, alignmentSet):
        if not alignmentSet.endswith(".alignmentset.xml"):
            raise InvalidDataset, "%s not an alignmentset"
        else:
            return alignmentSet
