#!/usr/bin/env python3
"""Offline snapshot contract fixtures, runnable with Python 3.8; no bindings."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

try:
    import h5py
    import numpy as np
except ImportError as error:
    print('SKIP: snapshot tool tests require numpy and h5py: '+str(error))
    sys.exit(77)

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from chui_snapshot import SnapshotError, read_snapshot, sha256, species_index, qualification_physics, normalized_hash
import hierarchical_dmo_qualification
import treepm_sweep


class SnapshotToolsTest(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory(prefix='chui_snapshot_tools_')
        self.root = Path(self.workspace.name)
        self.paths = []

    def tearDown(self):
        self.workspace.cleanup()

    def fixture(self, counts=(2, 1), mass_lane=False, dialect='chui_native'):
        """Construct the production v6 header and v2 manifest, not a checkpoint."""
        self.paths = []
        for i, count in enumerate(counts):
            path = self.root/('snap_0007.hdf5' if len(counts) == 1 else 'snap_0007.%d.hdf5' % i)
            self.paths.append(path)
            with h5py.File(path, 'w') as f:
                f.attrs['cosmosim_file_kind'] = 'science_snapshot'
                h = f.create_group('Header').attrs
                local = np.zeros(6, dtype='u8'); local[1] = count
                total = np.zeros(6, dtype='u8'); total[1] = sum(counts)
                for key, value in dict(NumPart_ThisFile=local, NumPart_Total=total,
                                      NumFilesPerSnapshot=np.uint32(len(counts)),
                                      CHUISnapshotMemberIndex=np.uint32(i), Time=0.5, Redshift=1.0,
                                      BoxSize=1.0 if dialect == 'chui_native' else 0.7,
                                      CHUIBoxSizeX_MpcComoving=1.0, CHUIBoxSizeY_MpcComoving=2.0,
                                      CHUIBoxSizeZ_MpcComoving=3.0, Omega0=0.3, OmegaLambda=0.7,
                                      OmegaBaryon=0.05, HubbleParam=0.7, CosmoSimSchemaVersion=np.uint32(6),
                                      CosmoSimSchemaName='chui_science_snapshot_v6',
                                      CHUISnapshotDialect=dialect, CHUISnapshotGenerationID='generation_7',
                                      NamingRulesVersion='2', FileNamingRulesVersion='2',
                                      CHUIVelocityStorageConvention=('chui_internal_peculiar_velocity' if
                                          dialect == 'chui_native' else 'peculiar_velocity_div_sqrt_a')).items():
                    h[key] = value
                h['MassTable'] = np.array([0, 2, 0, 0, 0, 0], dtype=float)
                units = f.create_group('Units').attrs
                for key, value in dict(LengthUnit='mpc', MassUnit='msun', VelocityUnit='km_s',
                                      CoordinateFrame='comoving', ExternalDialect=dialect).items():
                    units[key] = value
                text = 'schema_version = 1\n[cosmology]\nomega_matter = 0.3\n'
                provenance = f.create_group('Provenance').attrs
                provenance['config_hash_hex'] = normalized_hash(text)
                provenance['normalized_config_hash_hex'] = normalized_hash(text)
                provenance['normalized_config_sha256_hex'] = hashlib.sha256(text.encode()).hexdigest()
                provenance['normalized_config'] = text
                f.create_group('Config').attrs['normalized'] = text
                f.create_group('Parameters').attrs['normalized_param_txt'] = text
                block = f.create_group('PartType1')
                # Deliberately reverse local order and leave IDs noncontiguous.
                ids = np.arange(sum(counts[:i])+1, sum(counts[:i+1])+1, dtype='u8')[::-1]*11
                block['ParticleIDs'] = ids
                block['Coordinates'] = np.repeat((ids/100.)[:, None], 3, axis=1)
                block['Velocities'] = np.zeros((count, 3))
                if mass_lane:
                    block['Masses'] = np.full(count, 2.)
        return self.marker()

    def marker(self, version=2):
        with h5py.File(self.paths[0], 'r') as f:
            h = f['Header'].attrs
            fields = dict(schema='chui_snapshot_set_v%d' % version, generation_id=h['CHUISnapshotGenerationID'],
                          num_files_per_snapshot=int(h['NumFilesPerSnapshot']), schema_name=h['CosmoSimSchemaName'],
                          schema_version=int(h['CosmoSimSchemaVersion']), file_kind='science_snapshot',
                          dialect=h['CHUISnapshotDialect'], global_part_count=','.join(str(int(x)) for x in h['NumPart_Total']),
                          time=h['Time'], redshift=h['Redshift'], box_size_x=1, box_size_y=2, box_size_z=3,
                          omega_matter=h['Omega0'], omega_lambda=h['OmegaLambda'], omega_baryon=h['OmegaBaryon'],
                          hubble_param=h['HubbleParam'], unit_length='mpc', unit_mass='msun', unit_velocity='km_s',
                          coordinate_frame='comoving', velocity_storage_convention=h['CHUIVelocityStorageConvention'],
                          config_hash_hex=f['Provenance'].attrs['config_hash_hex'], naming_rules_version='2', file_naming_rules_version='2',
                          member_count=len(self.paths))
        for i, path in enumerate(self.paths):
            prefix = 'member.%d.' % i
            with h5py.File(path, 'r') as f:
                h = f['Header'].attrs
                fields.update({prefix+'filename': path.name, prefix+'index': int(h['CHUISnapshotMemberIndex']),
                               prefix+'local_part_count': ','.join(str(int(x)) for x in h['NumPart_ThisFile'])})
            fields[prefix+'file_size'] = path.stat().st_size
            fields[prefix+('sha256' if version == 2 else 'integrity')] = (
                sha256(path) if version == 2 else 'distributed_science_readback_v1')
        marker = self.root/'snap_0007.complete'
        body = ''.join('%s=%s\n' % item for item in fields.items())
        marker.write_text(body+'set_digest_sha256='+hashlib.sha256(body.encode()).hexdigest()+'\n')
        return marker

    def mutate(self, action, member=0):
        with h5py.File(self.paths[member], 'r+') as f:
            action(f)
        self.marker()

    def rejects(self, category, path=None):
        with self.assertRaisesRegex(SnapshotError, '^'+category+':'):
            read_snapshot(path or self.root/'snap_0007.complete', chunk_rows=1)

    def test_complete_set_mass_table_identity_and_both_consumers(self):
        marker = self.fixture()
        state = read_snapshot(self.paths[1], chunk_rows=1)
        np.testing.assert_array_equal(state['ids'], [11, 22, 33])
        np.testing.assert_array_equal(state['mass'], [2, 2, 2])
        np.testing.assert_array_equal(state['box'], [1, 2, 3])
        np.testing.assert_allclose(state['x'][:, 0], state['ids']/100., rtol=0, atol=0)
        self.assertEqual(state['provenance']['total_particle_count'], 3)
        self.assertEqual(len(state['provenance']['members']), 2)
        self.assertEqual(state['provenance']['manifest_sha256'], sha256(marker))
        self.assertEqual(hierarchical_dmo_qualification.read(marker)['provenance'], state['provenance'])
        class Args:
            snapshot = marker
            group = 'PartType1'
            legacy_dialect = legacy_length_unit = legacy_mass_unit = legacy_velocity_unit = None
            output = self.root/'particles.csv'
            target_ids = self.root/'targets.txt'
            epsilon = 0.001
            position_scale = mass_scale = 1.
            targets = 2
        treepm_sweep.export_snapshot(Args())
        self.assertEqual(Args.target_ids.read_text(), '11\n22\n')
        self.assertTrue(Args.output.with_suffix('.metadata.json').is_file())

    def test_fixed_targets_preserve_oracle_order_and_reject_bad_identity(self):
        marker = self.fixture()
        ids = read_snapshot(marker)['ids']
        targets = self.root/'historical_targets.txt'
        targets.write_text('33\n11\n')
        self.assertEqual(treepm_sweep.select_target_ids(ids, 2, targets), [33, 11])
        for text in ('11\n11\n', '44\n', '0\n', '-1\n', '11.0\n',
                     str(2**64)+'\n', '', '11\n\n'):
            targets.write_text(text)
            with self.assertRaisesRegex(ValueError, '^targets:'):
                treepm_sweep.select_target_ids(ids, 2, targets)

    def test_total_reference_requires_matching_source_and_independent_evidence(self):
        particles, targets, reference = (self.root/name for name in ('p.csv', 'ids.txt', 'r.csv'))
        particles.write_text('frozen source fixture\n')
        targets.write_text('11\n33\n')
        reference.write_text('id,total_x,total_y,total_z\n11,1,0,0\n33,-1,0,0\n')
        manifest = self.root/'reference.json'
        record = dict(reference_total_sha256=sha256(reference), particles_sha256=sha256(particles),
                      target_ids_sha256=sha256(targets), scale_factor=.5, box_solver_units=[1, 2, 3],
                      epsilon_comoving=.001, g_code=1., reference_kind='independent_periodic_ewald',
                      qualification_evidence='fixture: independent solver and convergence record')
        def check(value):
            manifest.write_text(json.dumps(value))
            return treepm_sweep.validate_total_reference(manifest, reference, particles, targets,
                                                        .5, [1, 2, 3], .001, 1.)
        self.assertEqual(check(record), record)
        for key, value in (('particles_sha256', 'wrong'), ('target_ids_sha256', 'wrong'),
                           ('scale_factor', .6), ('box_solver_units', [1, 1, 1]),
                           ('epsilon_comoving', .002), ('g_code', 2.),
                           ('reference_total_sha256', 'wrong'), ('reference_kind', 'same_solver'),
                           ('qualification_evidence', '')):
            bad = dict(record); bad[key] = value
            with self.assertRaisesRegex(ValueError, '^reference:'):
                check(bad)

    def test_single_file_v3_and_interoperable_geometry(self):
        self.fixture((3,), mass_lane=True, dialect='arepo_format3')
        marker = self.marker(version=3)
        state = read_snapshot(marker)
        np.testing.assert_allclose(state['box'], [.7, 1.4, 2.1], rtol=1e-15, atol=0)
        self.assertEqual(state['provenance']['identity']['stored_to_internal_position'], 1/.7)
        self.assertEqual(state['provenance']['identity']['stored_to_internal_velocity'], np.sqrt(.5))

    def test_missing_manifest_and_legacy_override_cannot_bypass_modern_closure(self):
        self.fixture((3,)).unlink()
        self.rejects('manifest', self.paths[0])
        with self.assertRaisesRegex(SnapshotError, '^manifest:'):
            read_snapshot(self.paths[0], legacy_dialect='chui_native')

    def test_missing_member_and_checksum_corruption(self):
        self.fixture(); self.paths[1].unlink()
        self.rejects('members')

    def test_bad_manifest_digest(self):
        marker = self.fixture()
        marker.write_text(marker.read_text().replace('member_count=2', 'member_count=3'))
        self.rejects('manifest')

    def test_member_hash_mismatch(self):
        self.fixture()
        with h5py.File(self.paths[1], 'r+') as f:
            f['PartType1/Coordinates'][0, 0] += 0.1
        self.rejects('manifest')

    def test_member_identity_conflicts(self):
        for key, value in [('Time', .6), ('CHUISnapshotGenerationID', 'other'),
                           ('CHUIBoxSizeY_MpcComoving', 4.), ('HubbleParam', .6)]:
            with self.subTest(key=key):
                self.fixture()
                self.mutate(lambda f: f['Header'].attrs.__setitem__(key, value), 1)
                self.rejects('identity')

    def test_duplicate_member_index(self):
        self.fixture()
        self.mutate(lambda f: f['Header'].attrs.__setitem__('CHUISnapshotMemberIndex', np.uint32(0)), 1)
        self.rejects('members')

    def test_global_count_mismatch(self):
        self.fixture()
        for i in range(2):
            self.mutate(lambda f: f['Header'].attrs.__setitem__('NumPart_Total', np.array([0, 4, 0, 0, 0, 0], dtype='u8')), i)
        self.rejects('counts')

    def test_duplicate_particle_ids_across_members(self):
        self.fixture()
        self.mutate(lambda f: f['PartType1/ParticleIDs'].__setitem__(0, 11), 1)
        self.rejects('ids')

    def test_group_contract(self):
        for group in ('PartType-1', 'PartType6', 'PartType01', 'PartType', 'Other1', '../PartType1'):
            with self.subTest(group=group), self.assertRaisesRegex(SnapshotError, '^group:'):
                species_index(group)
        self.assertEqual([species_index('PartType%d' % i) for i in range(6)], list(range(6)))

    def test_mass_and_population_rejections(self):
        for value in (None, [0]*6, [0, -1, 0, 0, 0, 0], [0, float('nan'), 0, 0, 0, 0], [1]):
            with self.subTest(table=value):
                self.fixture()
                def change(f):
                    del f['Header'].attrs['MassTable']
                    if value is not None:
                        f['Header'].attrs['MassTable'] = value
                self.mutate(change)
                self.rejects('mass')
        self.fixture((0,)); self.rejects('population')

    def test_invalid_lanes_units_and_checkpoint(self):
        for lane, value, category in [('Coordinates', float('nan'), 'lanes'),
                                      ('Velocities', float('inf'), 'lanes'), ('Masses', 0., 'lanes')]:
            with self.subTest(lane=lane):
                self.fixture(mass_lane=True)
                self.mutate(lambda f: f['PartType1/'+lane].__setitem__(0, value))
                self.rejects(category)
        self.fixture()
        self.mutate(lambda f: f['Units'].attrs.__setitem__('LengthUnit', 'parsec_guess'))
        self.rejects('units')
        self.fixture()
        self.mutate(lambda f: f.attrs.__setitem__('cosmosim_file_kind', 'restart_checkpoint'))
        self.rejects('file_kind')

    def test_unit_metadata_conflicts(self):
        for action in (
                lambda f: f['Header'].attrs.__setitem__('CHUIUnitLength', 'kpc'),
                lambda f: f['Header'].attrs.__setitem__('UnitLength_in_cm', 1.0),
                lambda f: f['Units'].attrs.__setitem__('CoordinateFrame', 'physical'),
                lambda f: f['Units'].attrs.__setitem__('ExternalDialect', 'gadget4_hdf5'),
                lambda f: f['PartType1/Coordinates'].attrs.__setitem__('h_scaling', -1.0)):
            self.fixture(); self.mutate(action); self.rejects('units')

    def test_invalid_id_lane_and_missing_velocity(self):
        self.fixture()
        self.mutate(lambda f: f['PartType1/ParticleIDs'].__setitem__(0, 0))
        self.rejects('ids')
        self.fixture()
        self.mutate(lambda f: f['PartType1'].__delitem__('Velocities'))
        self.rejects('lanes')

    def test_duplicate_manifest_member_path_and_order(self):
        for original, replacement, category in (
                ('member.1.filename=snap_0007.1.hdf5', 'member.1.filename=snap_0007.0.hdf5', 'members'),
                ('member.1.index=1', 'member.1.index=0', 'manifest')):
            marker = self.fixture()
            body = marker.read_text().split('set_digest_sha256=')[0].replace(original, replacement)
            marker.write_text(body+'set_digest_sha256='+hashlib.sha256(body.encode()).hexdigest()+'\n')
            self.rejects(category)


    def test_configuration_provenance_fingerprints_and_representations(self):
        for group, key, value in (('Provenance', 'config_hash_hex', 'wrong'),
                                  ('Provenance', 'normalized_config_hash_hex', 'wrong'),
                                  ('Provenance', 'normalized_config_sha256_hex', '0'*64),
                                  ('Provenance', 'normalized_config', 'different'),
                                  ('Config', 'normalized', 'different')):
            self.fixture((3,))
            self.mutate(lambda f: f[group].attrs.__setitem__(key, value))
            self.rejects('provenance')

    def test_scientific_configuration_comparison_is_constrained(self):
        text = ('schema_version = 1\n[cosmology]\nomega_matter = 0.3\n'
                '[numerics]\nhierarchical_max_rung = 0\ngravity_softening = 1 kpc\n'
                'treepm_pm_grid_nx = 64\n[output]\nrun_name = global\n')
        self.assertEqual(qualification_physics(text), qualification_physics(
            text.replace('hierarchical_max_rung = 0', 'hierarchical_max_rung = 2').replace('global', 'refined')))
        for old, new in (('omega_matter = 0.3', 'omega_matter = 0.4'),
                         ('gravity_softening = 1 kpc', 'gravity_softening = 2 kpc'),
                         ('treepm_pm_grid_nx = 64', 'treepm_pm_grid_nx = 128')):
            self.assertNotEqual(qualification_physics(text), qualification_physics(text.replace(old, new)))


    def test_explicit_external_legacy_path(self):
        self.fixture((3,)).unlink()
        with h5py.File(self.paths[0], 'r+') as f:
            del f.attrs['cosmosim_file_kind']
            for key in ('CosmoSimSchemaName', 'CosmoSimSchemaVersion', 'CHUISnapshotDialect'):
                del f['Header'].attrs[key]
        state = read_snapshot(self.paths[0], legacy_dialect='chui_native')
        self.assertEqual(state['provenance']['manifest_schema'], 'explicit_external_legacy')


if __name__ == '__main__':
    unittest.main()
