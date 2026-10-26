#!/usr/bin/env python3
"""Identity-translation regression for tools/compare_dmo_arepo.py (2026-10-26).

AUTHORED, NOT RUN. Pure-function tests (no snapshot files): explicit
--arepo-id-offset normalization, strict one-to-one matching and halo matching
on normalized identity. Exit 77 (ctest SKIP) when numpy/h5py are unavailable,
because the tool imports the validated snapshot reader.
"""
from pathlib import Path
import sys
import unittest

try:
    import numpy as np
    import h5py  # noqa: F401  (required by tools/chui_snapshot)
except ImportError as error:
    print('SKIP: compare_dmo_arepo identity tests require numpy and h5py: ' + str(error))
    sys.exit(77)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
import compare_dmo_arepo as tool  # noqa: E402

U64_MAX = (1 << 64) - 1


class IdentityNormalizationTest(unittest.TestCase):
    def test_identical_zero_based_sets(self):
        ids = np.arange(0, 8, dtype=np.uint64)
        oc, oa, matched, record = tool.match_identities(ids, ids.copy(), 0)
        np.testing.assert_array_equal(matched, ids)
        np.testing.assert_array_equal(ids[oc], ids[oa])
        self.assertEqual(record['arepo_id_offset'], 0)
        self.assertEqual(record['matched_particles'], 8)

    def test_chui_one_based_versus_arepo_zero_based(self):
        chui = np.arange(1, 65, dtype=np.uint64)
        arepo = np.arange(0, 64, dtype=np.uint64)
        oc, oa, matched, record = tool.match_identities(chui, arepo, 1)
        np.testing.assert_array_equal(chui[oc], arepo[oa] + np.uint64(1))
        self.assertEqual(record['normalized_id_range'], [1, 64])
        self.assertEqual(record['arepo_original_id_range'], [0, 63])
        self.assertIn('arepo_id_offset', record['convention'])
        # Originals are not modified.
        np.testing.assert_array_equal(arepo, np.arange(0, 64, dtype=np.uint64))

    def test_shuffled_order_matches_by_identity(self):
        rng = np.random.default_rng(7)
        base = np.arange(100, 132, dtype=np.uint64)
        chui = base + np.uint64(1)
        arepo = base.copy()
        rng.shuffle(chui)
        rng.shuffle(arepo)
        oc, oa, matched, _ = tool.match_identities(chui, arepo, 1)
        np.testing.assert_array_equal(chui[oc], arepo[oa] + np.uint64(1))
        np.testing.assert_array_equal(matched, np.sort(chui))

    def test_missing_particle_is_rejected(self):
        chui = np.arange(1, 11, dtype=np.uint64)
        arepo = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8], dtype=np.uint64)  # ID 9 missing
        with self.assertRaises(tool.IdentityMatchError) as ctx:
            tool.match_identities(chui, arepo, 1)
        self.assertIn('only_in_chui=1', str(ctx.exception))

    def test_duplicate_particle_is_rejected(self):
        chui = np.array([1, 2, 3, 3], dtype=np.uint64)
        arepo = np.array([0, 1, 2, 3], dtype=np.uint64)
        with self.assertRaises(tool.IdentityMatchError) as ctx:
            tool.match_identities(chui, arepo, 1)
        self.assertIn('duplicate', str(ctx.exception))
        with self.assertRaises(tool.IdentityMatchError):
            tool.match_identities(np.array([1, 2, 3, 4], dtype=np.uint64),
                                  np.array([0, 1, 1, 3], dtype=np.uint64), 1)

    def test_incorrect_offset_is_rejected_with_suggestion_but_not_applied(self):
        chui = np.arange(1, 11, dtype=np.uint64)
        arepo = np.arange(0, 10, dtype=np.uint64)
        with self.assertRaises(tool.IdentityMatchError) as ctx:
            tool.match_identities(chui, arepo, 0)  # default offset: no silent inference
        self.assertIn('--arepo-id-offset 1 would match', str(ctx.exception))
        with self.assertRaises(tool.IdentityMatchError):
            tool.match_identities(chui, arepo, 2)

    def test_very_large_64_bit_ids(self):
        top = np.uint64(U64_MAX)
        arepo = np.array([U64_MAX - 3, U64_MAX - 2, U64_MAX - 1], dtype=np.uint64)
        chui = np.array([U64_MAX - 2, U64_MAX - 1, U64_MAX], dtype=np.uint64)
        oc, oa, matched, record = tool.match_identities(chui, arepo, 1)
        self.assertEqual(int(matched[-1]), U64_MAX)
        self.assertEqual(record['normalized_id_range'][1], U64_MAX)
        self.assertEqual(int(top), U64_MAX)

    def test_overflow_and_underflow_are_rejected(self):
        with self.assertRaises(tool.IdentityMatchError) as ctx:
            tool.normalize_ids(np.array([U64_MAX], dtype=np.uint64), 1, 'arepo')
        self.assertIn('overflows', str(ctx.exception))
        with self.assertRaises(tool.IdentityMatchError) as ctx:
            tool.normalize_ids(np.array([0, 5], dtype=np.uint64), -1, 'arepo')
        self.assertIn('underflows', str(ctx.exception))
        with self.assertRaises(tool.IdentityMatchError):
            tool.normalize_ids(np.array([-1, 2], dtype=np.int64), 0, 'arepo')
        with self.assertRaises(tool.IdentityMatchError):
            tool.normalize_ids(np.array([1.0, 2.0]), 0, 'arepo')  # floating-point IDs
        with self.assertRaises(tool.IdentityMatchError):
            tool.normalize_ids(np.array([1, 2], dtype=np.uint64), 1.0, 'arepo')  # non-integer offset
        with self.assertRaises(tool.IdentityMatchError):
            tool.normalize_ids(np.array([1, 2], dtype=np.uint64), True, 'arepo')  # bool is not an offset
        # Signed 64-bit input is accepted when non-negative and converted exactly.
        np.testing.assert_array_equal(
            tool.normalize_ids(np.array([0, 2**62], dtype=np.int64), 1, 'arepo'),
            np.array([1, 2**62 + 1], dtype=np.uint64))

    def test_equal_count_different_identities_is_rejected(self):
        chui = np.array([1, 2, 3, 4], dtype=np.uint64)
        arepo = np.array([10, 20, 30, 40], dtype=np.uint64)
        with self.assertRaises(tool.IdentityMatchError) as ctx:
            tool.match_identities(chui, arepo, 1)
        self.assertNotIn('would match', str(ctx.exception))


class HaloMatchingTest(unittest.TestCase):
    def test_halo_membership_uses_normalized_identity(self):
        box = 1.0
        n = 40
        # Two clumps of 20 particles; AREPO zero-based and shuffled.
        chui_ids = np.arange(1, n + 1, dtype=np.uint64)
        arepo_ids = np.arange(0, n, dtype=np.uint64)
        positions = np.zeros((n, 3))
        positions[:20] = [0.25, 0.25, 0.25]
        positions[20:] = [0.75, 0.75, 0.75]
        positions += 1e-3 * np.arange(n)[:, None] / n
        rng = np.random.default_rng(3)
        perm = rng.permutation(n)
        arepo_ids_shuffled = arepo_ids[perm]
        arepo_positions = positions[perm]
        oc, oa, matched, _ = tool.match_identities(chui_ids, arepo_ids_shuffled, 1)
        xc = positions[oc]
        xa = arepo_positions[oa]
        np.testing.assert_allclose(xc, xa)  # same particle at the same matched index
        labels_c = np.where(np.arange(n) < 20, 0, 1)
        labels_a = labels_c.copy()
        mass = np.ones(n)
        hc = tool.halo_catalogue(labels_c, xc, mass, box, 5)
        ha = tool.halo_catalogue(labels_a, xa, mass, box, 5)
        matches = tool.match_halo_catalogues(hc, ha, matched, box, box / n ** (1.0 / 3.0))
        self.assertEqual(len(matches), 2)
        for match in matches:
            self.assertEqual(match['shared_members'], 20)
            self.assertAlmostEqual(match['mass_ratio'], 1.0)
            self.assertAlmostEqual(match['offset_over_spacing'], 0.0)


if __name__ == '__main__':
    unittest.main()
