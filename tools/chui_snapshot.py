"""Bounded offline reader for the production CHUI snapshot-set v2/v3 contract.

Python >=3.8. Science snapshots only: never a continuation checkpoint reader.
Arrays are returned in stored units; explicit caller conversions are recorded.
Native units and cosmological GADGET/AREPO storage conventions are validated,
not inferred from filenames. Legacy external inputs require an explicit dialect
and unit triple. No production configuration parser lives here.
"""
import hashlib
import math
from pathlib import Path
import re

import h5py
import numpy as np


class SnapshotError(ValueError):
    """A logical snapshot failed a named integrity/semantic contract."""


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def normalized_hash(text):
    """The production stable FNV-1a identity; SHA-256 remains separate."""
    value = 14695981039346656037
    for byte in text.encode('utf-8'):
        value = ((value ^ byte)*1099511628211) & ((1 << 64)-1)
    return format(value, '016x')


def species_index(group):
    if not isinstance(group, str) or re.fullmatch(r'PartType[0-5]', group) is None:
        raise SnapshotError('group: expected PartType0 through PartType5')
    return int(group[len('PartType'):])


def _string(value):
    if isinstance(value, bytes):
        return value.decode('utf-8').rstrip('\0')
    if isinstance(value, str):
        return value.rstrip('\0')
    raise SnapshotError('metadata: expected scalar string')


def _uint(value, label):
    array = np.asarray(value)
    if array.shape != () or array.dtype.kind not in 'ui' or int(array) < 0:
        raise SnapshotError('counts: invalid integer ' + label)
    return int(array)


def _counts(value, label):
    array = np.asarray(value)
    if array.shape != (6,) or array.dtype.kind not in 'ui' or (array < 0).any():
        raise SnapshotError('counts: expected six nonnegative integers for ' + label)
    return tuple(int(x) for x in array)


def _finite(value, label, positive=False):
    if np.asarray(value).shape != () or np.asarray(value).dtype.kind not in 'fiu':
        raise SnapshotError('metadata: expected scalar ' + label)
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0):
        raise SnapshotError('metadata: invalid finite ' + label)
    return result


def _equal(a, b):
    if isinstance(a, float) and isinstance(b, float):
        return abs(a-b) <= 64*np.finfo(float).eps*max(1, abs(a), abs(b))
    return a == b


def _manifest(path):
    if path.stat().st_size > 4*1024*1024:
        raise SnapshotError('manifest: exceeds 4 MiB metadata budget')
    raw = path.read_bytes()
    prefix, separator, tail = raw.rpartition(b'set_digest_sha256=')
    if not separator or re.fullmatch(b'[0-9a-f]{64}\n?', tail) is None:
        raise SnapshotError('manifest: missing/invalid set digest')
    if hashlib.sha256(prefix).hexdigest() != tail.decode().strip():
        raise SnapshotError('manifest: set digest mismatch')
    fields = {}
    for line in raw.decode('utf-8').splitlines():
        if not line:
            continue
        key, sep, value = line.partition('=')
        if not sep or not key or key in fields:
            raise SnapshotError('manifest: malformed or duplicate key')
        fields[key] = value
    if fields.get('schema') not in ('chui_snapshot_set_v2', 'chui_snapshot_set_v3'):
        raise SnapshotError('manifest: unsupported schema')
    return fields, hashlib.sha256(raw).hexdigest()


def _manifest_int(fields, key):
    value = fields.get(key, '')
    if re.fullmatch(r'[0-9]+', value) is None:
        raise SnapshotError('manifest: invalid integer ' + key)
    return int(value)


def _manifest_counts(fields, key):
    value = fields.get(key, '')
    if re.fullmatch(r'[0-9]+(?:,[0-9]+){5}', value) is None:
        raise SnapshotError('manifest: invalid count vector ' + key)
    return tuple(int(x) for x in value.split(','))


def _stem(path):
    return re.sub(r'(?:\.[0-9]+)?\.hdf5$', '', path.name)


def qualification_physics(text):
    """Compare canonical provenance; never interpret it as simulation input.

    Only output/run labels, execution placement and named time-resolution keys
    may differ between matched global/hierarchical convergence runs.
    """
    allowed_sections = ('output', 'parallel')
    allowed = {'run.run_name', 'run.output_root', 'run.log_level',
               'numerics.hierarchical_max_rung', 'numerics.max_global_steps',
               'numerics.cosmology_max_delta_ln_a', 'numerics.cosmology_max_hubble_time_fraction',
               'numerics.a_end', 'numerics.z_end', 'numerics.t_code_end', 'numerics.t_phys_end'}
    section = ''
    result = {}
    seen = set()
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith(('#', ';')):
            continue
        if line.startswith('[') and line.endswith(']'):
            section = line[1:-1]
            continue
        key, sep, value = line.partition('=')
        if not sep or (not section and key.strip() != 'schema_version'):
            raise SnapshotError('provenance: malformed normalized configuration')
        full = section+'.'+key.strip() if section else key.strip()
        if full in seen:
            raise SnapshotError('provenance: duplicate normalized configuration key')
        seen.add(full)
        if section not in allowed_sections and full not in allowed:
            result[full] = value.strip()
    if not result:
        raise SnapshotError('provenance: empty scientific configuration')
    return result


def _header(path, legacy_dialect, legacy_units):
    with h5py.File(path, 'r') as source:
        kind = _string(source.attrs.get('cosmosim_file_kind', ''))
        if kind not in ('', 'science_snapshot') or 'integrator' in source or 'state' in source:
            raise SnapshotError('file_kind: science snapshot required, checkpoint refused')
        if 'Header' not in source:
            raise SnapshotError('header: missing Header')
        h = source['Header'].attrs
        for key in ('NumPart_ThisFile', 'NumPart_Total', 'NumFilesPerSnapshot', 'Time', 'BoxSize'):
            if key not in h:
                raise SnapshotError('header: missing ' + key)
        local = _counts(h['NumPart_ThisFile'], 'NumPart_ThisFile')
        total = _counts(h['NumPart_Total'], 'NumPart_Total')
        if np.asarray(h['NumPart_Total']).dtype.itemsize == 4 and 'NumPart_Total_HighWord' in h:
            high = _counts(h['NumPart_Total_HighWord'], 'NumPart_Total_HighWord')
            if np.asarray(h['NumPart_Total_HighWord']).dtype.itemsize != 4:
                raise SnapshotError('counts: HighWord must be 32-bit')
            total = tuple(lo | (hi << 32) for lo, hi in zip(total, high))
        count = _uint(h['NumFilesPerSnapshot'], 'NumFilesPerSnapshot')
        if not 1 <= count <= 65536:
            raise SnapshotError('members: invalid count or member budget exceeded')
        match = re.fullmatch(r'(.+)\.([0-9]+)\.hdf5', path.name)
        index = _uint(h['CHUISnapshotMemberIndex'], 'member index') if 'CHUISnapshotMemberIndex' in h else (
            int(match.group(2)) if match else 0)
        if index >= count or (match and index != int(match.group(2))):
            raise SnapshotError('members: filename/header index disagreement')
        if count > 1 and not match and 'CHUISnapshotMemberIndex' not in h:
            raise SnapshotError('members: missing member index')
        epoch = _finite(h['Time'], 'Time', True)
        scalar_box = np.asarray(h['BoxSize'], dtype=float)
        if scalar_box.shape == ():
            box = np.repeat(scalar_box, 3)
        elif scalar_box.shape == (3,):
            box = scalar_box.copy()
        else:
            raise SnapshotError('geometry: unsupported BoxSize shape')
        axis_keys = ['CHUIBoxSize'+axis+'_MpcComoving' for axis in 'XYZ']
        if any(key in h for key in axis_keys):
            if not all(key in h for key in axis_keys) or scalar_box.shape != ():
                raise SnapshotError('geometry: incomplete/ambiguous CHUI axis geometry')
            physical_box = np.array([_finite(h[key], key, True) for key in axis_keys])
            box = float(scalar_box)*physical_box/physical_box[0]
        if not np.isfinite(box).all() or (box <= 0).any():
            raise SnapshotError('geometry: invalid box lengths')
        schema = _string(h.get('CosmoSimSchemaName', ''))
        version = _uint(h.get('CosmoSimSchemaVersion', np.uint32(0)), 'schema version')
        if schema.startswith('chui_science_snapshot') and kind != 'science_snapshot':
            raise SnapshotError('file_kind: modern CHUI science identity is incomplete')
        if kind == 'science_snapshot' and not (
                (schema == 'chui_science_snapshot_v6' and version == 6) or
                (schema == 'gadget_arepo_v'+str(version) and 1 <= version <= 5)):
            raise SnapshotError('schema: unsupported CHUI science snapshot revision')
        stored_dialect = _string(h.get('CHUISnapshotDialect', ''))
        if stored_dialect and legacy_dialect and stored_dialect != legacy_dialect:
            raise SnapshotError('dialect: legacy override conflicts with stored dialect')
        dialect = stored_dialect or legacy_dialect
        if dialect not in ('chui_native', 'arepo_format3', 'gadget4_hdf5'):
            raise SnapshotError('dialect: explicit supported legacy dialect required')
        if not stored_dialect and kind == 'science_snapshot' and (not schema or version > 5):
            raise SnapshotError('dialect: malformed modern CHUI snapshot')
        units_group = source['Units'].attrs if 'Units' in source else {}
        if 'ExternalDialect' in units_group and _string(units_group['ExternalDialect']) != dialect:
            raise SnapshotError('units: ExternalDialect disagrees with snapshot dialect')
        unit_keys = ('LengthUnit', 'MassUnit', 'VelocityUnit')
        header_keys = ('CHUIUnitLength', 'CHUIUnitMass', 'CHUIUnitVelocity')
        units = []
        for i, (key, header_key) in enumerate(zip(unit_keys, header_keys)):
            value = units_group.get(key, h.get(header_key, ''))
            value = _string(value)
            if key in units_group and header_key in h and value != _string(h[header_key]):
                raise SnapshotError('units: Header/Units disagreement')
            if value and legacy_units and value != legacy_units[i]:
                raise SnapshotError('units: explicit legacy units conflict with stored units')
            units.append(value or (legacy_units[i] if legacy_units else ''))
        for value, supported in zip(units, (('m', 'cm', 'kpc', 'mpc'), ('kg', 'g', 'msun'), ('m_s', 'cm_s', 'km_s'))):
            if value not in supported:
                raise SnapshotError('units: missing or unsupported unit convention')
        cgs_units = (dict(m=100.0, cm=1.0, kpc=3.0856775814913673e21, mpc=3.0856775814913673e24)[units[0]],
                     dict(kg=1000.0, g=1.0, msun=1.98847e33)[units[1]],
                     dict(m_s=100.0, cm_s=1.0, km_s=1e5)[units[2]])
        for attributes in (h, units_group):
            for key, expected_cgs in zip(('UnitLength_in_cm', 'UnitMass_in_g', 'UnitVelocity_in_cm_per_s'), cgs_units):
                if key in attributes and not _equal(_finite(attributes[key], key, True), expected_cgs):
                    raise SnapshotError('units: CGS factors disagree with declared unit convention')
        frame = _string(units_group.get('CoordinateFrame', 'comoving' if legacy_dialect else ''))
        if frame != 'comoving':
            raise SnapshotError('units: these cosmological tools require comoving coordinates')
        convention = _string(h.get('CHUIVelocityStorageConvention', '')) or (
            'chui_internal_peculiar_velocity' if dialect == 'chui_native' else 'peculiar_velocity_div_sqrt_a')
        expected = 'chui_internal_peculiar_velocity' if dialect == 'chui_native' else 'peculiar_velocity_div_sqrt_a'
        if convention != expected:
            raise SnapshotError('units: unsupported velocity storage convention')
        provenance = source['Provenance'].attrs if 'Provenance' in source else {}
        identity = dict(kind=kind, schema=schema, version=version, dialect=dialect, units=tuple(units),
                        coordinate_frame=frame, velocity_storage_convention=convention,
                        generation=_string(h.get('CHUISnapshotGenerationID', '')),
                        config_hash=_string(provenance.get('config_hash_hex', '')),
                        normalized_config_hash=_string(provenance.get('normalized_config_hash_hex', '')),
                        normalized_config_sha256=_string(provenance.get('normalized_config_sha256_hex', '')),
                        naming_rules=_string(h.get('NamingRulesVersion', '')),
                        file_naming_rules=_string(h.get('FileNamingRulesVersion', '')),
                        time=epoch, box=tuple(float(x) for x in box), total=total, count=count)
        if identity['generation'] and (Path(identity['generation']).name != identity['generation'] or '\\' in identity['generation'] or identity['generation'] in ('.', '..')):
            raise SnapshotError('identity: unsafe generation identifier')
        for key in ('Redshift', 'Omega0', 'OmegaLambda', 'OmegaBaryon', 'HubbleParam'):
            identity[key] = _finite(h[key], key) if key in h else None
        if dialect != 'chui_native' and (identity['HubbleParam'] is None or identity['HubbleParam'] <= 0):
            raise SnapshotError('units: positive HubbleParam required for cosmological storage')
        for key in ('CHUIGravitySourceGeneration', 'CHUISourceGenerationID', 'CHUISourceIdentity'):
            if key in h:
                value = h[key]
                identity[key] = _string(value) if isinstance(value, (str, bytes)) else _uint(value, key)
        if 'normalized_param_txt' in (source['Parameters'].attrs if 'Parameters' in source else {}):
            text = _string(source['Parameters'].attrs['normalized_param_txt'])
            if len(text.encode('utf-8')) > 4*1024*1024:
                raise SnapshotError('provenance: normalized config exceeds metadata budget')
            identity['normalized_text_sha256'] = hashlib.sha256(text.encode()).hexdigest()
            expected_hash = normalized_hash(text)
            for key in ('config_hash', 'normalized_config_hash'):
                if identity[key] and identity[key] != expected_hash:
                    raise SnapshotError('provenance: normalized configuration fingerprint mismatch')
            if identity['normalized_config_sha256'] and identity['normalized_config_sha256'] != identity['normalized_text_sha256']:
                raise SnapshotError('provenance: normalized configuration SHA-256 mismatch')
            for attributes, key in ((provenance, 'normalized_config'),
                    (source['Config'].attrs if 'Config' in source else {}, 'normalized')):
                if key in attributes and _string(attributes[key]) != text:
                    raise SnapshotError('provenance: normalized configuration representations disagree')
            identity['qualification_physics'] = qualification_physics(text)
        if all(key in h for key in axis_keys):
            length_si = dict(m=1.0, cm=0.01, kpc=3.0856775814913673e19,
                             mpc=3.0856775814913673e22)[units[0]]
            stored_per_mpc = 3.0856775814913673e22/length_si
            if dialect != 'chui_native':
                stored_per_mpc *= identity['HubbleParam']
            expected_box = physical_box*stored_per_mpc
            if not np.allclose(box, expected_box, rtol=64*np.finfo(float).eps, atol=0):
                raise SnapshotError('geometry: stored BoxSize disagrees with physical axes/units/dialect')
        identity['stored_to_internal_position'] = 1.0 if dialect == 'chui_native' else 1.0/identity['HubbleParam']
        identity['stored_to_internal_mass'] = identity['stored_to_internal_position']
        identity['stored_to_internal_velocity'] = 1.0 if dialect == 'chui_native' else math.sqrt(epoch)
    return dict(path=path, local=local, index=index, identity=identity)


def _validate_marker(marker, fields, members):
    ref = members[0]['identity']
    if marker.stem not in (ref['generation'], _stem(members[0]['path'])):
        raise SnapshotError('manifest: marker stem disagrees with generation/logical snapshot')
    if _manifest_int(fields, 'member_count') != len(members) or _manifest_int(fields, 'num_files_per_snapshot') != ref['count']:
        raise SnapshotError('manifest: member count mismatch')
    if _manifest_counts(fields, 'global_part_count') != ref['total']:
        raise SnapshotError('manifest: global counts mismatch')
    strings = {'generation_id': 'generation', 'schema_name': 'schema', 'file_kind': 'kind',
               'dialect': 'dialect', 'coordinate_frame': 'coordinate_frame', 'config_hash_hex': 'config_hash',
               'velocity_storage_convention': 'velocity_storage_convention',
               'naming_rules_version': 'naming_rules', 'file_naming_rules_version': 'file_naming_rules'}
    for key, identity_key in strings.items():
        if key not in fields or fields[key] != ref[identity_key]:
            raise SnapshotError('manifest: identity mismatch ' + key)
    if _manifest_int(fields, 'schema_version') != ref['version']:
        raise SnapshotError('manifest: schema version mismatch')
    for key, value in zip(('unit_length', 'unit_mass', 'unit_velocity'), ref['units']):
        if fields.get(key) != value:
            raise SnapshotError('manifest: units mismatch ' + key)
    numbers = dict(time=ref['time'], redshift=ref['Redshift'], omega_matter=ref['Omega0'],
                   omega_lambda=ref['OmegaLambda'], omega_baryon=ref['OmegaBaryon'], hubble_param=ref['HubbleParam'])
    # Production manifests store physical Mpc axis extensions, rather than the
    # interoperable BoxSize unit. Read those exact values for marker validation.
    with h5py.File(members[0]['path'], 'r') as source:
        h = source['Header'].attrs
        for axis in 'xyz':
            numbers['box_size_'+axis] = _finite(h.get('CHUIBoxSize'+axis.upper()+'_MpcComoving', ref['box']['xyz'.index(axis)]), 'box')
    for key, expected in numbers.items():
        try:
            actual = float(fields[key])
        except (KeyError, ValueError):
            raise SnapshotError('manifest: invalid scientific field ' + key)
        if expected is None or not math.isfinite(actual) or not _equal(actual, expected):
            raise SnapshotError('manifest: scientific mismatch ' + key)
    if fields['schema'] == 'chui_snapshot_set_v3' and len(members) != 1:
        raise SnapshotError('manifest: v3 distributed readback requires one file')
    for i, member in enumerate(members):
        prefix = 'member.'+str(i)+'.'
        path = member['path']
        if (fields.get(prefix+'filename') != path.name or _manifest_int(fields, prefix+'index') != member['index'] or
                _manifest_counts(fields, prefix+'local_part_count') != member['local'] or
                _manifest_int(fields, prefix+'file_size') != path.stat().st_size):
            raise SnapshotError('manifest: member identity/size/count mismatch')
        if fields['schema'] == 'chui_snapshot_set_v2':
            if fields.get(prefix+'sha256') != sha256(path):
                raise SnapshotError('manifest: member SHA-256 mismatch')
        elif fields.get(prefix+'integrity') != 'distributed_science_readback_v1':
            raise SnapshotError('manifest: missing distributed science readback evidence')


def _validate_lanes(source, local, identity):
    table = source['Header'].attrs.get('MassTable')
    if table is not None:
        if np.asarray(table).dtype.kind not in 'fiu':
            raise SnapshotError('mass: nonnumeric MassTable')
        table = np.asarray(table, dtype=float)
        if table.shape != (6,) or not np.isfinite(table).all() or (table < 0).any():
            raise SnapshotError('mass: invalid MassTable')
    for part, n in enumerate(local):
        key = 'PartType'+str(part)
        if key not in source:
            if n:
                raise SnapshotError('lanes: missing populated ' + key)
            continue
        block = source[key]
        for lane, shape in (('ParticleIDs', (n,)), ('Coordinates', (n, 3)), ('Velocities', (n, 3))):
            if lane not in block or block[lane].shape != shape:
                raise SnapshotError('lanes: missing/invalid shape ' + key+'/'+lane)
        interop = identity['dialect'] != 'chui_native'
        cgs = (dict(m=100.0, cm=1.0, kpc=3.0856775814913673e21, mpc=3.0856775814913673e24)[identity['units'][0]],
               dict(m_s=100.0, cm_s=1.0, km_s=1e5)[identity['units'][2]],
               dict(kg=1000.0, g=1.0, msun=1.98847e33)[identity['units'][1]])
        for lane, a_power, h_power, factor in zip(('Coordinates', 'Velocities', 'Masses'),
                (1.0, 0.5 if interop else 0.0, 0.0), (-1.0 if interop else 0.0, 0.0, -1.0 if interop else 0.0), cgs):
            if lane not in block:
                continue
            if block[lane].dtype.kind not in 'fiu':
                raise SnapshotError('lanes: nonnumeric '+key+'/'+lane)
            for attr, expected in (('a_scaling', a_power), ('h_scaling', h_power), ('CGSConversionFactor', factor)):
                if attr in block[lane].attrs and not _equal(_finite(block[lane].attrs[attr], attr), expected):
                    raise SnapshotError('units: dataset scaling conflicts with declared dialect/units')
        if block['ParticleIDs'].dtype.kind not in 'ui' or block['ParticleIDs'].dtype.itemsize not in (4, 8):
            raise SnapshotError('ids: expected 32/64-bit integer ParticleIDs')
        if 'Masses' in block and block['Masses'].shape != (n,):
            raise SnapshotError('mass: invalid Masses shape')
        if n and 'Masses' not in block and (table is None or table[part] <= 0):
            raise SnapshotError('mass: missing Masses or positive MassTable entry')
    return table


def read_snapshot(path, group='PartType1', legacy_dialect=None, legacy_units=None,
                  need_velocity=True, chunk_rows=65536):
    """Read a complete logical set, validating all particle groups in chunks.

    Only selected species XYZ/mass/(optional velocity) arrays are retained.
    One ID lane across all species is necessary to reject cross-member/type
    identity collisions. Sorting reorders selected arrays one lane at a time.
    """
    species = species_index(group)
    if chunk_rows <= 0:
        raise SnapshotError('budget: chunk_rows must be positive')
    path = Path(path).resolve()
    if not path.exists():
        raise SnapshotError('manifest: completion marker missing' if path.suffix == '.complete' else 'members: input missing')
    marker = None
    fields = None
    marker_checksum = None
    if path.is_dir():
        markers = sorted(path.glob('*.complete'))
        if len(markers) != 1:
            raise SnapshotError('manifest: directory requires exactly one completion marker')
        marker = markers[0]
    elif path.suffix == '.complete':
        marker = path
    else:
        first = _header(path, legacy_dialect, legacy_units)
        candidates = [path.parent/(_stem(path)+'.complete')]
        if first['identity']['generation']:
            candidates.append(path.parent/(first['identity']['generation']+'.complete'))
        marker = next((p for p in candidates if p.is_file()), None)
    if marker is not None:
        if not marker.is_file():
            raise SnapshotError('manifest: completion marker missing')
        fields, marker_checksum = _manifest(marker)
        count = _manifest_int(fields, 'member_count')
        if not 1 <= count <= 65536:
            raise SnapshotError('members: invalid manifest member budget')
        indices = {int(re.match(r'member\.([0-9]+)\.', key).group(1)) for key in fields if re.match(r'member\.([0-9]+)\.', key)}
        if indices != set(range(count)):
            raise SnapshotError('manifest: missing or extraneous member index')
        paths = []
        for i in range(count):
            name = fields.get('member.'+str(i)+'.filename', '')
            if not name or Path(name).name != name or '\\' in name or name in ('.', '..'):
                raise SnapshotError('manifest: unsafe member filename')
            member = marker.parent/name
            if not member.is_file():
                raise SnapshotError('members: missing member ' + name)
            paths.append(member)
        if len(set(p.resolve() for p in paths)) != len(paths):
            raise SnapshotError('members: duplicate member path')
        if path.is_file() and path.suffix != '.complete' and path not in [p.resolve() for p in paths]:
            raise SnapshotError('manifest: requested member absent from marker')
    else:
        if not legacy_dialect:
            raise SnapshotError('manifest: missing .complete; external legacy dialect must be explicit')
        if first['identity']['kind'] == 'science_snapshot' and not (
                re.fullmatch(r'gadget_arepo_v[1-5]', first['identity']['schema']) is not None and 1 <= first['identity']['version'] <= 5):
            raise SnapshotError('manifest: modern CHUI science snapshot requires .complete')
        paths = [path] if first['identity']['count'] == 1 else sorted(path.parent.glob(_stem(path)+'.*.hdf5'))
    initial_stats = {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in paths}
    members = [_header(p, legacy_dialect, legacy_units) for p in paths]
    ref = members[0]['identity']
    if len(members) != ref['count']:
        raise SnapshotError('members: incomplete NumFilesPerSnapshot population')
    if set(m['index'] for m in members) != set(range(ref['count'])):
        raise SnapshotError('members: duplicate or missing member index')
    if fields is not None and [m['index'] for m in members] != list(range(ref['count'])):
        raise SnapshotError('members: manifest ordering disagrees with member indices')
    members.sort(key=lambda m: m['index'])
    if len(set(_stem(m['path']) for m in members)) != 1:
        raise SnapshotError('members: conflicting logical snapshot stems')
    for member in members:
        current = member['identity']
        if current.keys() != ref.keys() or any(not _equal(current[k], ref[k]) for k in ref):
            raise SnapshotError('identity: conflicting epoch/generation/geometry/units/provenance')
    if tuple(sum(m['local'][i] for m in members) for i in range(6)) != ref['total']:
        raise SnapshotError('counts: member totals disagree with global totals')
    if marker is not None:
        _validate_marker(marker, fields, members)
    population = ref['total'][species]
    if population == 0:
        raise SnapshotError('population: selected particle species is empty')
    # Check on-disk extents and dtypes before admitting any population array.
    for member in members:
        with h5py.File(member['path'], 'r') as source:
            _validate_lanes(source, member['local'], member['identity'])
    all_ids = np.empty(sum(ref['total']), dtype=np.uint64)
    ids = np.empty(population, dtype=np.uint64)
    xyz = np.empty((population, 3), dtype=float)
    mass = np.empty(population, dtype=float)
    velocity = np.empty((population, 3), dtype=float) if need_velocity else None
    selected_offset = 0
    identity_offset = 0
    for member in members:
        with h5py.File(member['path'], 'r') as source:
            table = _validate_lanes(source, member['local'], member['identity'])
            for part in range(6):
                n = member['local'][part]
                key = 'PartType'+str(part)
                if key not in source:
                    continue
                block = source[key]
                for begin in range(0, n, chunk_rows):
                    end = min(begin+chunk_rows, n)
                    raw_ids = block['ParticleIDs'][begin:end]
                    if (raw_ids < 0).any() or (ref['kind'] == 'science_snapshot' and (raw_ids == 0).any()):
                        raise SnapshotError('ids: invalid persistent ParticleID')
                    pos = np.asarray(block['Coordinates'][begin:end], dtype=float)
                    vel = np.asarray(block['Velocities'][begin:end], dtype=float)
                    masses = (np.asarray(block['Masses'][begin:end], dtype=float) if 'Masses' in block else
                              np.full(end-begin, table[part]))
                    if not np.isfinite(pos).all() or not np.isfinite(vel).all() or not np.isfinite(masses).all() or (masses <= 0).any():
                        raise SnapshotError('lanes: nonfinite coordinates/velocities or nonpositive mass')
                    all_ids[identity_offset:identity_offset+end-begin] = raw_ids
                    identity_offset += end-begin
                    if part == species:
                        target = slice(selected_offset, selected_offset+end-begin)
                        ids[target] = raw_ids
                        xyz[target] = pos
                        mass[target] = masses
                        if velocity is not None:
                            velocity[target] = vel
                        selected_offset += end-begin
    all_ids.sort()
    if (all_ids[1:] == all_ids[:-1]).any():
        raise SnapshotError('ids: duplicate ParticleIDs across snapshot members/species')
    global_fingerprint = hashlib.sha256()
    for begin in range(0, len(all_ids), chunk_rows):
        global_fingerprint.update(all_ids[begin:begin+chunk_rows].astype('<u8', copy=False).tobytes())
    del all_ids
    order = np.argsort(ids, kind='stable')
    ids = ids[order]
    xyz = xyz[order]
    mass = mass[order]
    if velocity is not None:
        velocity = velocity[order]
    fingerprint = hashlib.sha256()
    for begin in range(0, population, chunk_rows):
        fingerprint.update(ids[begin:begin+chunk_rows].astype('<u8', copy=False).tobytes())
    provenance = dict(manifest_path=str(marker) if marker else None,
                      manifest_sha256=marker_checksum,
                      manifest_schema=fields['schema'] if fields else 'explicit_external_legacy',
                      members=[dict(path=str(m['path']), sha256=sha256(m['path']), index=m['index'],
                                    local_counts=m['local']) for m in members],
                      total_particle_count=sum(ref['total']), selected_particle_count=population,
                      particle_id_sha256_le_u64=fingerprint.hexdigest(),
                      all_species_id_sha256_le_u64=global_fingerprint.hexdigest(), identity=ref,
                      legacy_override=dict(dialect=legacy_dialect, units=legacy_units))
    if any((p.stat().st_size, p.stat().st_mtime_ns) != initial_stats[p] for p in paths):
        raise SnapshotError('identity: snapshot member changed while being consumed')
    if marker is not None and sha256(marker) != marker_checksum:
        raise SnapshotError('identity: completion marker changed while being consumed')
    return dict(ids=ids, x=xyz, u=velocity, mass=mass, box=np.array(ref['box']), a=ref['time'], provenance=provenance)


def add_snapshot_arguments(parser):
    parser.add_argument('--legacy-dialect', choices=('chui_native', 'arepo_format3', 'gadget4_hdf5'),
                        help='explicit external/CHUI <=v5 compatibility path; modern CHUI still requires .complete')
    parser.add_argument('--legacy-length-unit', choices=('m', 'cm', 'kpc', 'mpc'))
    parser.add_argument('--legacy-mass-unit', choices=('kg', 'g', 'msun'))
    parser.add_argument('--legacy-velocity-unit', choices=('m_s', 'cm_s', 'km_s'))


def legacy_units(args):
    values = (args.legacy_length_unit, args.legacy_mass_unit, args.legacy_velocity_unit)
    if any(values) and not all(values):
        raise SnapshotError('units: legacy overrides require the complete unit triple')
    return values if all(values) else None
