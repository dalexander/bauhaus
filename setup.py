from __future__ import absolute_import
from __future__ import print_function
from setuptools import setup, Extension, find_packages
import os.path
import sys

REQUIREMENTS_TXT = "requirements.txt"

globals = {}
exec(compile(open("bauhaus/__version__.py").read(), "bauhaus/__version__.py", 'exec'), globals)
__VERSION__ = globals["__VERSION__"]

def _get_local_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)

def _get_requirements(file_name):
    with open(file_name, 'r') as f:
        reqs = [line for line in f if not line.startswith("#")]
    return reqs

def _get_local_requirements(file_name):
    return _get_requirements(_get_local_file(file_name))

setup(
    name = "bauhaus",
    version=__VERSION__,
    author="David Alexander",
    author_email="dalexander@pacificbiosciences.com",
    description="A simplistic tertiary analysis engine",
    license=open("COPYING").read(),
    packages = find_packages('.'),
    zip_safe = False,
    entry_points = {
        "console_scripts" : [ "bauhaus = bauhaus.main:main" ]
    },
    package_data={ "bauhaus.resources" : [ "run.sh",
                                           "conditions.json",
                                           "tool-contracts/*.json",
                                           "scripts/R/*.R" ] },
    install_requires=_get_local_requirements(REQUIREMENTS_TXT)
)
