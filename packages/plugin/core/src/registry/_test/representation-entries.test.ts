/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { BuiltInStructureRepresentations } from '@molstar/graphics/repr/structure/catalog';
import { BuiltInVolumeRepresentations } from '@molstar/graphics/repr/volume/catalog';
import { BuiltInParticleRepresentations } from '@molstar/graphics/repr/particles/catalog';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { setProductionMode } from '@molstar/core/util/debug';
import { PluginContext } from '@molstar/plugin/context';
import {
  DefaultParticleRepresentations,
  DefaultStructureRepresentations,
  DefaultVolumeRepresentations,
} from '@molstar/plugin/default-registry';
import { Backbone } from '@molstar/plugin/registry/structure/backbone';
import { BallAndStick } from '@molstar/plugin/registry/structure/ball-and-stick';
import { BlobSurface } from '@molstar/plugin/registry/structure/blob-surface';
import { Carbohydrate } from '@molstar/plugin/registry/structure/carbohydrate';
import { Cartoon } from '@molstar/plugin/registry/structure/cartoon';
import { Ellipsoid } from '@molstar/plugin/registry/structure/ellipsoid';
import { GaussianSurface } from '@molstar/plugin/registry/structure/gaussian-surface';
import { GaussianVolume } from '@molstar/plugin/registry/structure/gaussian-volume';
import { Label } from '@molstar/plugin/registry/structure/label';
import { Line } from '@molstar/plugin/registry/structure/line';
import { MolecularSurface } from '@molstar/plugin/registry/structure/molecular-surface';
import { Orientation } from '@molstar/plugin/registry/structure/orientation';
import { Plane } from '@molstar/plugin/registry/structure/plane';
import { Point } from '@molstar/plugin/registry/structure/point';
import { Polyhedron } from '@molstar/plugin/registry/structure/polyhedron';
import { Putty } from '@molstar/plugin/registry/structure/putty';
import { Spacefill } from '@molstar/plugin/registry/structure/spacefill';
import { DirectVolume } from '@molstar/plugin/registry/volume/direct-volume';
import { Dot } from '@molstar/plugin/registry/volume/dot';
import { Isosurface } from '@molstar/plugin/registry/volume/isosurface';
import { Segment } from '@molstar/plugin/registry/volume/segment';
import { Slice } from '@molstar/plugin/registry/volume/slice';
import { ParticleFibers } from '@molstar/plugin/registry/particles/fibers';
import { ParticleOrientation } from '@molstar/plugin/registry/particles/orientation';
import { ParticleSpacefill } from '@molstar/plugin/registry/particles/spacefill';
import { ParticleTarget } from '@molstar/plugin/registry/particles/target';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

type Scope = 'structure' | 'volume' | 'particles';

const StructureEntries = [
  Cartoon,
  Backbone,
  BallAndStick,
  BlobSurface,
  Carbohydrate,
  Ellipsoid,
  GaussianSurface,
  GaussianVolume,
  Label,
  Line,
  MolecularSurface,
  Orientation,
  Plane,
  Point,
  Putty,
  Spacefill,
  Polyhedron,
];
const VolumeEntries = [DirectVolume, Dot, Isosurface, Segment, Slice];
const ParticleEntries = [ParticleSpacefill, ParticleOrientation, ParticleFibers, ParticleTarget];

const Entries: [Scope, PluginRegistryEntry[], readonly { name: string }[]][] = [
  ['structure', StructureEntries, Object.values(BuiltInStructureRepresentations)],
  ['volume', VolumeEntries, Object.values(BuiltInVolumeRepresentations)],
  ['particles', ParticleEntries, Object.values(BuiltInParticleRepresentations)],
];

describe('representation registry entries', () => {
  for (const [scope, entries, catalog] of Entries) {
    describe(scope, () => {
      it('has one entry per built-in representation, in catalog order', () => {
        const providers = entries.map((e) => e[scope]!.representations![0]);
        expect(providers).toEqual(catalog);
        providers.forEach((p, i) => expect(p).toBe(catalog[i]));
      });

      it("carries the provider's default color and size themes and nothing else", () => {
        for (const entry of entries) {
          expect(Object.keys(entry)).toEqual([scope]);
          const e = entry[scope]!;
          expect(Object.keys(e).sort()).toEqual(['representations', 'themes']);
          expect(e.representations!.length).toBe(1);
          const provider = e.representations![0];

          const { color, size } = e.themes!;
          expect(color!.map((t) => t.name)).toEqual([provider.defaultColorTheme.name]);
          expect(size!.map((t) => t.name)).toEqual([provider.defaultSizeTheme.name]);
          // The very providers of the built-in themes, not copies.
          expect(color![0]).toBe((BuiltInColorThemes as Record<string, unknown>)[provider.defaultColorTheme.name]);
          expect(size![0]).toBe((BuiltInSizeThemes as Record<string, unknown>)[provider.defaultSizeTheme.name]);
        }
      });
    });
  }

  it('the default representation entries list the same providers in the same order', () => {
    expect(DefaultStructureRepresentations.structure!.representations).toEqual(
      Object.values(BuiltInStructureRepresentations),
    );
    expect(DefaultVolumeRepresentations.volume!.representations).toEqual(Object.values(BuiltInVolumeRepresentations));
    expect(DefaultParticleRepresentations.particles!.representations).toEqual(
      Object.values(BuiltInParticleRepresentations),
    );
  });

  it('entry modules do not import catalogs', () => {
    const root = path.resolve(__dirname, '..');
    let count = 0;
    for (const scope of ['structure', 'volume', 'particles']) {
      for (const file of fs.readdirSync(path.join(root, scope))) {
        const source = fs.readFileSync(path.join(root, scope, file), 'utf8');
        const imports = source.split('\n').filter((l) => l.startsWith('import '));
        expect(imports.filter((l) => /catalog/.test(l))).toEqual([]);
        count++;
      }
    }
    expect(count).toBe(StructureEntries.length + VolumeEntries.length + ParticleEntries.length);
  });
});

describe('default theme check', () => {
  // the representation alone, without the entry that brings its themes
  const missing = [
    "Volume representation 'direct-volume' default color theme 'volume-value' is not registered",
    "Volume representation 'direct-volume' default size theme 'uniform' is not registered",
  ];

  async function createPlugin(registry: PluginRegistryEntry[] = []) {
    const plugin = new PluginContext({ behaviors: [], registry });
    const warnings: string[] = [];
    plugin.events.log.subscribe((e) => {
      if (e.type === 'warning') warnings.push(e.message);
    });
    await plugin.init();
    return { plugin, warnings };
  }

  afterEach(() => setProductionMode(false));

  it('is silent for complete entries', async () => {
    const { warnings } = await createPlugin([...StructureEntries, ...VolumeEntries, ...ParticleEntries]);
    expect(warnings).toEqual([]);
  });

  it('is silent for the default plugin', async () => {
    const plugin = new PluginContext(DefaultPluginSpec());
    const warnings: string[] = [];
    plugin.events.log.subscribe((e) => {
      if (e.type === 'warning') warnings.push(e.message);
    });
    await plugin.init();
    expect(warnings).toEqual([]);
  });

  it('warns after register for a representation whose default theme is not registered, once', async () => {
    const { plugin, warnings } = await createPlugin();
    expect(warnings).toEqual([]);

    plugin.register({ volume: { representations: DirectVolume.volume!.representations } });
    expect(warnings).toEqual(missing);

    // Not repeated by later checks
    plugin.register(Isosurface);
    expect(warnings).toEqual(missing);

    // The check does not register anything
    expect(plugin.representation.volume.themes.colorThemeRegistry.has('volume-value')).toBe(false);
  });

  it('warns at the end of init for the spec registry', async () => {
    const plugin = new PluginContext({
      behaviors: [],
      registry: [{ volume: { representations: DirectVolume.volume!.representations } }],
    });
    const warnings: string[] = [];
    plugin.events.log.subscribe((e) => {
      if (e.type === 'warning') warnings.push(e.message);
    });
    await plugin.init();
    expect(warnings).toEqual(missing);
  });

  it('is silent when the entry brings the default themes', async () => {
    const { plugin, warnings } = await createPlugin();
    plugin.register(DirectVolume);
    expect(warnings).toEqual([]);
  });

  it('does not run in production mode', async () => {
    const { plugin, warnings } = await createPlugin();
    setProductionMode(true);
    plugin.register({ volume: { representations: DirectVolume.volume!.representations } });
    expect(warnings).toEqual([]);
  });
});
