from __future__ import annotations

import pathlib
import tempfile
import unittest

import cosmosim


class SnapshotDiagnosticIntegrationTest(unittest.TestCase):
    def test_roundtrip_snapshot_and_health(self) -> None:
        config = cosmosim.SimulationConfig()
        config.run_name = "python_bindings_test"

        with tempfile.TemporaryDirectory(prefix="cosmosim_py_") as tmp_dir:
            snapshot_path = pathlib.Path(tmp_dir) / "snapshot_python_000.hdf5"
            state = cosmosim.make_uniform_dark_matter_state(32, 1.0, config.box_size_mpc_comov)
            cosmosim.write_snapshot(
                snapshot_path,
                state,
                config,
                normalized_config_text="schema_version = 1\n",
                git_sha="testsha",
            )

            loaded = cosmosim.read_snapshot(snapshot_path, config)
            self.assertEqual(loaded.state.particle_count, 32)

            summary = cosmosim.summarize_run_health(loaded, config)
            self.assertTrue(summary["ownership_invariants_ok"])
            self.assertTrue(summary["unique_particle_ids_ok"])
            self.assertEqual(summary["particle_count"], 32)
            self.assertAlmostEqual(cosmosim.snapshot_mass_sum(loaded), 32.0)


if __name__ == "__main__":
    unittest.main()
