from __future__ import annotations

import unittest

import cosmosim


class BindingsSmokeTest(unittest.TestCase):
    def test_import_and_state_views(self) -> None:
        config = cosmosim.SimulationConfig()
        self.assertEqual(config.schema_version, 1)

        state = cosmosim.make_uniform_dark_matter_state(8, 2.5, 10.0)
        self.assertEqual(state.particle_count, 8)

        ids = state.particle_id()
        self.assertEqual(ids.shape, (8,))
        self.assertFalse(ids.flags.writeable)
        self.assertEqual(int(ids[0]), 1)

        mass = state.mass_code()
        self.assertAlmostEqual(float(mass.sum()), 20.0)


if __name__ == "__main__":
    unittest.main()
