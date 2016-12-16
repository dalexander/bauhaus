from pkg_resources import Requirement, resource_filename
import os.path as op

def getScriptPath(fname):
    path = resource_filename(Requirement.parse("bauhaus"), op.join("bauhaus/scripts/", fname))
    if not op.exists(path):
        raise ValueError, "Invalid script: %s" % path
    else:
        return path
