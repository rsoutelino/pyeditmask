import unittest
import matplotlib.pyplot as plt

from pyeditmask import WW3Grid

class TestWW3Grid(unittest.TestCase):

    def test_read_meta(self):
        filename = '/source/ww3/tools/gridgen/output/nz_coastal.meta'
        ww3grid = WW3Grid(filename)
        print ww3grid.maskr.shape
        print ww3grid.latr.shape
        print ww3grid.lonr.shape


if __name__ == "__main__":
    unittest.main() 
