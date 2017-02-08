
from builtins import object
class Workflow(object):
    """
    Abstract base class for workflows.
    (I tried making this work with abc.ABCMeta, it didn't like the static methods)
    """
    @staticmethod
    def conditionTableType():
        raise NotImplementedError

    @staticmethod
    def name():
        raise NotImplementedError

    def generate(self, pflow, ct):
       raise NotImplementedError
