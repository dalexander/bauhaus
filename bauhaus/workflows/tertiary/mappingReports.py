__all__ = [ "MappingReportsWorkflow"]

import os.path as op

from bauhaus import Workflow
from bauhaus.experiment import (InputType, ResequencingConditionTable)
#from .datasetOps import *

def generateConditionsJSON(pflow, ct):
    """
    Generate a pbcommandR-compatible "condition table" JSON
    """
    pass

def generateToolContract(pflow, toolName):
    """
    Generate a tool contract for running one of the mapping report
    plotting tools from "internaltools"
    """
    pass

class MappingReportsWorkflow(Workflow):
    """
    Basic mapping---not chunked, just dead simple.
    """
    @staticmethod
    def name():
        return "MappingReports"

    @staticmethod
    def conditionTableType():
        return ResequencingConditionTable

    def generate(self, pflow, ct):
        outputDict = {}
        for condition in ct.conditions:
            with pflow.context("condition", condition):
                pass
