import os
from unittest import TestCase

from .. import run

class testRun(TestCase):
    def test_config_parse(self):
        self.assertEqual(
	    run.config_parse(os.path.join("tests", "config_test.cfg")),
	    {"BLAST":
                {"blast":
                     "C:\\Program files\\NCBI\\blast-2.2.28+\\bin",
		 "blastp":
                     "C:\\Program files\\NCBI\\blast-2.2.28+\\bin\\blastp.exe"}
		})
