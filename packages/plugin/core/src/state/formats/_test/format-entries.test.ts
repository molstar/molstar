/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultActions, DefaultFormats } from '@molstar/plugin/default-registry';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { Isosurface } from '@molstar/plugin/registry/volume/isosurface';
import { Segment } from '@molstar/plugin/registry/volume/segment';
import { ParticleSpacefill } from '@molstar/plugin/registry/particles/spacefill';
import { ParticleFibers } from '@molstar/plugin/registry/particles/fibers';
import { ParticleTarget } from '@molstar/plugin/registry/particles/target';
import { ParseCif } from '../cif.js';
import { BuiltInVolumeFormats } from '../volume/catalog.js';
import { BuiltInTopologyFormats } from '../topology/catalog.js';
import { BuiltInCoordinatesFormats } from '../coordinates/catalog.js';
import { BuiltInShapeFormats } from '../shape/catalog.js';
import { BuiltInParticlesFormats } from '../particles/catalog.js';
import { BuiltInTrajectoryFormats } from '../trajectory/catalog.js';
import { Ccp4, ParseCcp4, VolumeFromCcp4 } from '../volume/ccp4.js';
import { Dsn6, ParseDsn6, VolumeFromDsn6 } from '../volume/dsn6.js';
import { Cube, VolumeFromCube } from '../volume/cube.js';
import { Dx, VolumeFromDx } from '../volume/dx.js';
import { Dscif } from '../volume/density-server.js';
import { Segcif } from '../volume/segmentation.js';
import { Sfcif } from '../volume/structure-factors.js';
import { Mtz } from '../volume/mtz.js';
import { Psf } from '../topology/psf.js';
import { Prmtop } from '../topology/prmtop.js';
import { Top } from '../topology/top.js';
import { Dcd } from '../coordinates/dcd.js';
import { Xtc } from '../coordinates/xtc.js';
import { Trr } from '../coordinates/trr.js';
import { Nctraj } from '../coordinates/nctraj.js';
import { LammpsTrajectory } from '../coordinates/lammps.js';
import { Ply } from '../shape/ply.js';
import { Obj } from '../shape/obj.js';
import { Vtp } from '../shape/vtp.js';
import { ParticleListFromRelionStar, RelionStarParticles } from '../particles/star.js';
import { DynamoTblParticles, ParticleListFromDynamoTbl } from '../particles/tbl.js';
import { CryoEtDataPortalNdjsonParticles, ParticleListFromCryoEtDataPortalNdjson } from '../particles/ndjson.js';
import { ArtiatomiEmParticles, ParticleListFromArtiatomiEm } from '../particles/em.js';
import { MmcifParticles, ParticleListFromMmcifAssembly } from '../particles/mmcif-assembly.js';
import { ParticleTrajectoryFromSimularium, SimulariumParticles } from '../particles/simularium.js';
import { Mmcif, TrajectoryFromMmCif } from '../trajectory/mmcif.js';
import { CifCore, TrajectoryFromCifCore } from '../trajectory/cif-core.js';
import { Pdb, Pdbqt, Pqr, TrajectoryFromPDB } from '../trajectory/pdb.js';
import { Gro } from '../trajectory/gro.js';
import { Xyz } from '../trajectory/xyz.js';
import { LammpsData, LammpsTrajectoryData } from '../trajectory/lammps.js';
import { Mol } from '../trajectory/mol.js';
import { Sdf } from '../trajectory/sdf.js';
import { Mol2 } from '../trajectory/mol2.js';
import { ComplexParticleVisuals, SimpleParticleVisuals } from '../particles/provider.js';

type Spec = {
  entry: PluginRegistryEntry;
  /** The built-in provider the entry registers. */
  provider: { name: string };
  /** The actions `DefaultActions` lists for the format, in order. */
  actions?: unknown[];
  /** The representation entry the format's visuals use. */
  visuals?: PluginRegistryEntry;
};

// In the order of `DefaultFormats`: volume, topology, coordinates, shape, particles, trajectory.
const Volume: Spec[] = [
  { entry: Ccp4, provider: BuiltInVolumeFormats[0], actions: [ParseCcp4, VolumeFromCcp4], visuals: Isosurface },
  { entry: Dsn6, provider: BuiltInVolumeFormats[1], actions: [ParseDsn6, VolumeFromDsn6], visuals: Isosurface },
  { entry: Cube, provider: BuiltInVolumeFormats[2], actions: [VolumeFromCube], visuals: Isosurface },
  { entry: Dx, provider: BuiltInVolumeFormats[3], actions: [VolumeFromDx], visuals: Isosurface },
  { entry: Dscif, provider: BuiltInVolumeFormats[4], visuals: Isosurface },
  { entry: Segcif, provider: BuiltInVolumeFormats[5], visuals: Segment },
  { entry: Sfcif, provider: BuiltInVolumeFormats[6], visuals: Isosurface },
  { entry: Mtz, provider: BuiltInVolumeFormats[7], visuals: Isosurface },
];
const Topology: Spec[] = [
  { entry: Psf, provider: BuiltInTopologyFormats[0] },
  { entry: Prmtop, provider: BuiltInTopologyFormats[1] },
  { entry: Top, provider: BuiltInTopologyFormats[2] },
];
const Coordinates: Spec[] = [
  { entry: Dcd, provider: BuiltInCoordinatesFormats[0] },
  { entry: Xtc, provider: BuiltInCoordinatesFormats[1] },
  { entry: Trr, provider: BuiltInCoordinatesFormats[2] },
  { entry: Nctraj, provider: BuiltInCoordinatesFormats[3] },
  { entry: LammpsTrajectory, provider: BuiltInCoordinatesFormats[4] },
];
const Shape: Spec[] = [
  { entry: Ply, provider: BuiltInShapeFormats[0] },
  { entry: Obj, provider: BuiltInShapeFormats[1] },
  { entry: Vtp, provider: BuiltInShapeFormats[2] },
];
const Particles: Spec[] = [
  {
    entry: RelionStarParticles,
    provider: BuiltInParticlesFormats[0],
    actions: [ParticleListFromRelionStar],
    visuals: SimpleParticleVisuals,
  },
  {
    entry: DynamoTblParticles,
    provider: BuiltInParticlesFormats[1],
    actions: [ParticleListFromDynamoTbl],
    visuals: SimpleParticleVisuals,
  },
  {
    entry: CryoEtDataPortalNdjsonParticles,
    provider: BuiltInParticlesFormats[2],
    actions: [ParticleListFromCryoEtDataPortalNdjson],
    visuals: SimpleParticleVisuals,
  },
  {
    entry: ArtiatomiEmParticles,
    provider: BuiltInParticlesFormats[3],
    actions: [ParticleListFromArtiatomiEm],
    visuals: SimpleParticleVisuals,
  },
  {
    entry: MmcifParticles,
    provider: BuiltInParticlesFormats[4],
    actions: [ParticleListFromMmcifAssembly],
    visuals: ComplexParticleVisuals,
  },
  {
    entry: SimulariumParticles,
    provider: BuiltInParticlesFormats[5],
    actions: [ParticleTrajectoryFromSimularium],
    visuals: ComplexParticleVisuals,
  },
];
const Trajectory: Spec[] = [
  { entry: Mmcif, provider: BuiltInTrajectoryFormats[0], actions: [ParseCif, TrajectoryFromMmCif] },
  { entry: CifCore, provider: BuiltInTrajectoryFormats[1], actions: [TrajectoryFromCifCore] },
  { entry: Pdb, provider: BuiltInTrajectoryFormats[2], actions: [TrajectoryFromPDB] },
  { entry: Pdbqt, provider: BuiltInTrajectoryFormats[3], actions: [TrajectoryFromPDB] },
  { entry: Pqr, provider: BuiltInTrajectoryFormats[4], actions: [TrajectoryFromPDB] },
  { entry: Gro, provider: BuiltInTrajectoryFormats[5] },
  { entry: Xyz, provider: BuiltInTrajectoryFormats[6] },
  { entry: LammpsData, provider: BuiltInTrajectoryFormats[7] },
  { entry: LammpsTrajectoryData, provider: BuiltInTrajectoryFormats[8] },
  { entry: Mol, provider: BuiltInTrajectoryFormats[9] },
  { entry: Sdf, provider: BuiltInTrajectoryFormats[10] },
  { entry: Mol2, provider: BuiltInTrajectoryFormats[11] },
];
const All = [...Volume, ...Topology, ...Coordinates, ...Shape, ...Particles, ...Trajectory];

describe('format entries', () => {
  it('has one entry per built-in format provider, in the order of DefaultFormats', () => {
    expect(All.length).toBe(DefaultFormats.formats!.length);
    All.forEach((spec, i) => {
      expect(spec.entry.formats).toEqual([DefaultFormats.formats![i]]);
      expect(spec.entry.formats![0]).toBe(spec.provider);
    });
  });

  it('carries the actions DefaultActions lists for the format and no others', () => {
    const defaults = new Set<unknown>(DefaultActions.actions!);
    for (const { entry, actions = [] } of All) {
      expect(entry.actions ?? []).toEqual(actions);
      // 5.x listed all of them, so DefaultRegistry registers no action that it did not
      for (const a of entry.actions ?? []) expect(defaults.has(a)).toBe(true);
    }
  });

  it('includes the representations and themes the visuals use, only for volume and particle formats', () => {
    for (const { entry, visuals } of All) {
      const { formats, actions, ...rest } = entry;
      if (!visuals) {
        // Topology, coordinates, shape, and trajectory formats carry no representation providers
        expect(rest).toEqual({});
        continue;
      }
      const { formats: _f, actions: _a, ...expected } = visuals;
      expect(rest).toEqual(expected);
    }
  });

  it('volume entries include the isosurface (segment for segmentations) and its uniform themes', () => {
    for (const { entry, visuals } of Volume) {
      expect(entry.volume).toBe(visuals!.volume);
    }
    const iso = Ccp4.volume!;
    expect(iso.representations!.map((r) => r.name)).toEqual(['isosurface']);
    expect(iso.themes!.color!.map((t) => t.name)).toEqual(['uniform']);
    expect(iso.themes!.size!.map((t) => t.name)).toEqual(['uniform']);
    expect(Segcif.volume!.representations!.map((r) => r.name)).toEqual(['segment']);
    expect(Segcif.volume!.themes!.color!.map((t) => t.name)).toEqual(['volume-segment']);
  });

  it('particle entries include what their visuals apply', () => {
    expect(RelionStarParticles.particles).toBe(ParticleSpacefill.particles);
    const complex = SimulariumParticles.particles!;
    expect(complex.representations!.map((r) => r.name)).toEqual(['spacefill', 'fibers', 'target']);
    expect(complex.themes!.color!.map((t) => t.name).sort()).toEqual([
      'particle-compartment',
      'particle-entity',
      'particle-hierarchy',
      'particle-index',
      'uniform',
    ]);
    // every default representation of the three particle entries is in the complex entry
    for (const e of [ParticleSpacefill, ParticleFibers, ParticleTarget]) {
      for (const r of e.particles!.representations!) expect(complex.representations).toContain(r);
      for (const t of e.particles!.themes!.color!) expect(complex.themes!.color).toContain(t);
      for (const t of e.particles!.themes!.size!) expect(complex.themes!.size).toContain(t);
    }
    expect(new Set(complex.themes!.size).size).toBe(complex.themes!.size!.length);
  });

  it('registers into a plugin', async () => {
    const plugin = new PluginContext({ behaviors: [], registry: [] });
    await plugin.init();
    expect(plugin.dataFormats.list).toEqual([]);
    const unregister = plugin.register([Sdf, Ccp4, Mmcif]);
    expect(plugin.dataFormats.has('sdf')).toBe(true);
    expect(plugin.dataFormats.has('ccp4')).toBe(true);
    expect(plugin.dataFormats.has('mmcif')).toBe(true);
    expect(plugin.dataFormats.list.map((f) => f.name)).toEqual(['sdf', 'ccp4', 'mmcif']);
    unregister();
    expect(plugin.dataFormats.list).toEqual([]);
    plugin.dispose();
  });

  it('format modules do not import catalogs or presets', () => {
    const root = path.resolve(__dirname, '..');
    let count = 0;
    for (const family of ['volume', 'topology', 'coordinates', 'shape', 'particles', 'trajectory']) {
      for (const file of fs.readdirSync(path.join(root, family))) {
        if (!file.endsWith('.ts') || file === 'catalog.ts') continue;
        const source = fs.readFileSync(path.join(root, family, file), 'utf8');
        const imports = source.split('\n').filter((l) => l.startsWith('import ') && !l.startsWith('import type'));
        expect(imports.filter((l) => /catalog|presets/.test(l))).toEqual([]);
        count++;
      }
    }
    expect(count).toBeGreaterThan(40);
  });
});
