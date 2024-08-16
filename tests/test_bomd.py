import unittest
from src.bomd import BOMD

class TestBOMD(unittest.TestCase):
    def setUp(self):
        self.bomd = BOMD("config.yaml")

    def test_initialization(self):
        self.bomd.initialize()
        self.assertIsNotNone(self.bomd.positions)
        self.assertIsNotNone(self.bomd.velocities)
        self.assertIsNotNone(self.bomd.atoms_list)
        self.assertEqual(self.bomd.timesteps, [0.0])

    # Add more test methods as needed

if __name__ == '__main__':
    unittest.main()