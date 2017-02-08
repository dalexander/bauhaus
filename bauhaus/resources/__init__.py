from pkg_resources import Requirement, resource_filename
import os.path as op

def getResourcePath(fname):
    path = resource_filename(Requirement.parse("bauhaus"), op.join("bauhaus/resources/", fname))
    if not op.exists(path):
        raise ValueError("Invalid resource: %s" % path)
    else:
        return path
