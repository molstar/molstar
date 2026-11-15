/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import { Asset } from '@molstar/core/util/assets';
import { CustomProperties } from '@molstar/model/model/custom-property';
import type { ParticleList, ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { writeParticlePack } from '@molstar/model/formats/particles/pack';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultRegistry } from '@molstar/plugin/default-registry';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { OpenFiles } from '../../../actions/file.js';
import { ParticlesRepresentation3D } from '../representation.js';
import { ParticleTileSetFromPack, ParticleTilesView } from '../tiles.js';

/** 200 particles over a box, two targets of sites alternating, an entity and a compartment each */
function createList(): ParticleList {
  const n = 200;
  const coordinates = new Float32Array(n * 3);
  for (let i = 0; i < n; ++i) coordinates.set([(i * 37) % 101, (i * 53) % 97, (i * 11) % 89], i * 3);
  const targets = Int32Array.from({ length: n }, (_, i) => (i % 2 ? 2 : 1));
  const sites = (count: number): ParticleTarget => ({
    kind: 'sites',
    sites: {
      count,
      positions: Float32Array.from({ length: count * 3 }, (_, k) => k % 3),
      radii: new Float32Array(count).fill(2),
    },
  });
  return {
    label: 'streamed',
    count: n,
    targets,
    targetInfo: new Map([
      [1, { name: 'one', entity: 0 }],
      [2, { name: 'two', entity: 1 }],
    ]),
    targetMapping: new Map([
      [1, sites(3)],
      [2, sites(2)],
    ]),
    entities: Int32Array.from(targets, (t) => t - 1),
    entityInfo: new Map([
      [0, { name: 'e0' }],
      [1, { name: 'e1' }],
    ]),
    compartments: Int32Array.from(targets, () => 0),
    compartmentInfo: new Map([[0, { name: 'c0' }]]),
    coordinates,
    rotations: new Float32Array(n * 4).map((_, k) => (k % 4 === 3 ? 1 : 0)),
    getParticleLabel: (i) => `p${i}`,
    sourceData: { kind: 'test', name: 'test', data: undefined },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

describe('a particle pack opened by range', () => {
  let plugin: PluginContext;
  beforeEach(async () => {
    plugin = new PluginContext({ behaviors: [], registry: [...DefaultRegistry, ...ParticlesRegistry] });
    await plugin.init();
  });
  afterEach(() => plugin.dispose());

  const open = (bytes: Uint8Array, name: string, visuals: boolean) =>
    plugin.runTask(
      plugin.state.data.applyAction(OpenFiles, {
        files: [Asset.File(new File([bytes], name))],
        format: { name: 'auto', params: {} },
        visuals,
      }),
    );
  const cellsOf = (transformer: unknown) => {
    const out: any[] = [];
    plugin.state.data.cells.forEach((cell) => {
      if (cell.transform.transformer === transformer) out.push(cell);
    });
    return out;
  };

  it('opens a tile-ordered pack as a tile set read by range, a view following the camera and its targets drawn', async () => {
    await open(writeParticlePack(createList(), { tiles: 4 }), 'cell.mspk', true);
    const [tiles] = cellsOf(ParticleTileSetFromPack);
    expect(tiles).toBeDefined();
    expect(SO.Particle.TileSet.is(tiles.obj)).toBe(true);
    expect(tiles.obj.data.set.count).toBe(200);
    expect(tiles.obj.description).toContain('read by range');
    // no list read whole
    let lists = 0;
    plugin.state.data.cells.forEach((cell) => {
      if (SO.Particle.List.is(cell.obj) && cell.transform.transformer !== ParticleTilesView) lists++;
    });
    expect(lists).toBe(0);
    const [view] = cellsOf(ParticleTilesView);
    expect(view.transform.params.tiles.name).toBe('camera');
    const [repr] = cellsOf(ParticlesRepresentation3D);
    expect(repr.transform.params.type.name).toBe('target');
    // entities and compartments: colored as a whole pack's
    expect(repr.transform.params.colorTheme.name).toBe('particle-hierarchy');
  });

  it('reads a pack by target whole, as before', async () => {
    await open(writeParticlePack(createList()), 'cell.mspk', false);
    expect(cellsOf(ParticleTileSetFromPack)).toHaveLength(0);
    let lists = 0;
    plugin.state.data.cells.forEach((cell) => {
      if (SO.Particle.List.is(cell.obj)) lists++;
    });
    expect(lists).toBe(1);
  });
});
