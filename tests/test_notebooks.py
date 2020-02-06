import os
import numpy as np
import testipynb
import unittest

NBDIR = os.path.sep.join(
    os.path.abspath(__file__).split(os.path.sep)[:-2] + ['notebooks']
)

IGNORE = ["2-power-line-plot.ipynb", "1-steel-casing-effect.ipynb"]

n_ignore = 5  # so we don't run over-time on travis, randomly ignore 5 notebooks
Test = testipynb.TestNotebooks(directory=NBDIR, timeout=2800)
ignore_inds = np.random.choice(len(Test._nbnames) - len(IGNORE), n_ignore)
test_nbnames = [t for t in Test._nbnames if t not in IGNORE]
Test.ignore = IGNORE + [test_nbnames[i] for i in ignore_inds]
TestNotebooks = Test.get_tests()

if __name__ == "__main__":
    unittest.main()
