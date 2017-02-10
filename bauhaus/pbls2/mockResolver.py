from builtins import object
__all__ = [ "MockResolver" ]

import os.path as op
from .resolver import _isRuncode
from .exceptions import *

class MockResolver(object):
    # For testing purposes

    REFERENCE_MASKS_ROOT = "/pbi/dept/consensus/bauhaus/genome-masks"
    REFERENCES_ROOT = "/mnt/secondary/iSmrtanalysis/current/common/references"

    def __init__(self):
        pass

    def resolveSubreadSet(self, runCode, reportsFolder=""):
        if not _isRuncode(runCode):
            raise ValueError('Argument "%s" does not appear to be a runcode' % runCode)
        lookup = \
            { ("3150128-0001", "")        : "/pbi/collections/315/3150128/r54008_20160308_001811/1_A01/m54008_160308_002050.subreadset.xml" ,
              ("3150128-0002", "")        : "/pbi/collections/315/3150128/r54008_20160308_001811/2_B01/m54008_160308_053311.subreadset.xml" ,
              ("3150128-0001", "Params1") : "/pbi/collections/315/3150128/r54008_20160308_001811/1_A01/Params1/m54008_160308_002050.subreadset.xml" ,
              ("3150128-0002", "Params1") : "/pbi/collections/315/3150128/r54008_20160308_001811/2_B01/Params1/m54008_160308_053311.subreadset.xml" ,
              ("3150128-0001", "Params2") : "/pbi/collections/315/3150128/r54008_20160308_001811/1_A01/Params2/m54008_160308_002050.subreadset.xml" ,
              ("3150128-0002", "Params2") : "/pbi/collections/315/3150128/r54008_20160308_001811/2_B01/Params2/m54008_160308_053311.subreadset.xml" ,
              ("3150122-0001", "")        : "/pbi/collections/315/3150122/r54011_20160305_235615/1_A01/m54011_160305_235923.subreadset.xml" ,
              ("3150122-0002", "")        : "/pbi/collections/315/3150122/r54011_20160305_235615/2_B01/m54011_160306_050740.subreadset.xml" }
        if (runCode, reportsFolder) not in lookup:
            raise DataNotFound("Input data not found: %s/%s" % (runCode, reportsFolder))
        return lookup[(runCode, reportsFolder)]

    def resolveReference(self, referenceName):
        if referenceName not in ["lambdaNEB", "ecoliK12_pbi_March2013", "plasmidbell_v1", "pBR322_EcoRV"]:
            raise DataNotFound("Reference not found: %s" % referenceName)
        referenceFasta = op.join(self.REFERENCES_ROOT, referenceName, "sequence", referenceName + ".fasta")
        return referenceFasta

    def resolveReferenceMask(self, referenceName):
        if referenceName not in ["lambdaNEB", "ecoliK12_pbi_March2013"]:
            raise DataNotFound("Reference mask (required for CoverageTitration) not found for genome: %s" % referenceName)
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
        if not (subreadSet.endswith(".subreadset.xml") or subreadSet.endswith(".subreads.bam")):
            raise InvalidDataset("%s not a subreadset" % subreadSet)
        else:
            return subreadSet

    def ensureAlignmentSet(self, alignmentSet):
        if not alignmentSet.endswith(".alignmentset.xml"):
            raise InvalidDataset("%s not an alignmentset")
        else:
            return alignmentSet

    def ensureTraceH5File(self, traceH5File):
        if not traceH5File.endswith(".trc.h5"):
            raise InvalidDataset("%s not an trc.h5 file" % traceH5File)
        else:
            return traceH5File
