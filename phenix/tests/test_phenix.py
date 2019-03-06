
from pyworkflow.tests import *
from phenix import Plugin
from phenix.constants import PHENIXVERSION

class TestVersion(BaseTest):

    def testgetVersion(self):
        version = Plugin.getPhenixVersion()
        self.assertEqual(version, PHENIXVERSION)
