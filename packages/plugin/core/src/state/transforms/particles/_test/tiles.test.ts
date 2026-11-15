/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Subject } from 'rxjs';
import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { OrderedSet } from '@molstar/core/data/int';
import {
  Particle,
  type ParticleList,
  type ParticleTarget,
  ParticleValues,
} from '@molstar/model/model/particles/particle-list';
import { ParticleTileView } from '@molstar/model/model/particles/particle-tile-view';
import { ParticleTiles } from '@molstar/model/model/particles/particle-tiles';
import { ArrayParticleTrajectory } from '@molstar/model/model/particles/particle-trajectory';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultRegistry } from '@molstar/plugin/default-registry';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { ParticleIndexColorThemeProvider } from '@molstar/graphics/theme/color/particle-index';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import {
  ParticleTileSetFromList,
  ParticleTileSetFromPack,
  ParticleTilesView,
  pickParticleTiles,
  setParticleTilesStep,
} from '@molstar/plugin/state/transforms/particles/tiles';
import {
  getParticleTilesStepParticles,
  ParticleTilesFollowCamera,
} from '@molstar/plugin/behavior/dynamic/particle-tiles';
import { PluginSpec } from '@molstar/plugin/spec';
import { StateSelection } from '@molstar/core/state';
import { Asset } from '@molstar/core/util/assets';

import { writeParticlePack } from '@molstar/model/formats/particles/pack';
import {
  createParticleTileSetFromPack,
  openParticlePackStream,
  ParticlePackSource,
} from '@molstar/model/formats/particles/pack-stream';

const Count = 64;

/** 64 particles on a 4³ lattice over 0-90 Å, targets 1 and 2 alternating along x (sites each) */
function createList(shift = 0): ParticleList {
  const coordinates = new Float32Array(Count * 3);
  const targets = new Int32Array(Count);
  for (let i = 0; i < Count; ++i) {
    const x = i % 4;
    coordinates.set([x * 30 + shift, (Math.floor(i / 4) % 4) * 30, Math.floor(i / 16) * 30], i * 3);
    targets[i] = x % 2 ? 2 : 1;
  }
  const sites = (n: number): ParticleTarget => ({
    kind: 'sites',
    sites: {
      count: n,
      positions: Float32Array.from({ length: n * 3 }, (_, k) => k % 3),
      radii: new Float32Array(n).fill(2),
    },
  });
  return {
    label: 'lattice',
    count: Count,
    targets,
    targetInfo: new Map([
      [1, { name: 'one' }],
      [2, { name: 'two' }],
    ]),
    targetMapping: new Map([
      [1, sites(2)],
      [2, sites(3)],
    ]),
    coordinates,
    getParticleLabel: (i) => `particle ${i}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

/** a list or trajectory put in the state as it is */
const ListOf = PluginStateTransform.BuiltIn({
  name: 'test-particle-list-of',
  from: SO.Root,
  to: SO.Particle.List,
  params: { list: PD.Value<ParticleList>(undefined as any, { isHidden: true }) },
})({
  apply: ({ params }) => new SO.Particle.List(params.list, { label: 'lattice' }),
});
const TrajectoryOf = PluginStateTransform.BuiltIn({
  name: 'test-particle-trajectory-of',
  from: SO.Root,
  to: SO.Particle.Trajectory,
  params: { trajectory: PD.Value<ArrayParticleTrajectory>(undefined as any, { isHidden: true }) },
})({
  apply: ({ params }) => new SO.Particle.Trajectory(params.trajectory, { label: 'lattice' }),
});

describe('particle tiles in the state', () => {
  let plugin: PluginContext;
  beforeEach(async () => {
    plugin = new PluginContext({ behaviors: [], registry: [...DefaultRegistry, ...ParticlesRegistry] });
    await plugin.init();
  });
  afterEach(() => plugin.dispose());

  async function tileSet(params: Partial<{ cells: number; delay: number; opacity: boolean }> = {}) {
    const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
    return plugin
      .build()
      .to(list.ref)
      .apply(ParticleTileSetFromList, { cells: 2, ...params })
      .commit();
  }

  it('subdivides a list into a tile set, a view of some tiles a particle list', async () => {
    const tiles = await tileSet();
    const { set, loader } = tiles.data!;
    expect(set.tiles.size).toBe(8);
    expect(set.count).toBe(Count);
    expect(tiles.cell!.obj!.description).toBe(`2³ grid · 8 tiles · ${Count} particles`);

    const view = await plugin
      .build()
      .to(tiles.ref)
      .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids: '0, 1' } }, frame: 0, slots: false })
      .commit();
    const list = view.data!;
    expect(list.count).toBe(set.tiles.get(0)!.count + set.tiles.get(1)!.count);
    expect(ParticleTileView.of(list)!.tiles).toEqual(Int32Array.of(0, 1));
    expect(loader.getState(0, 0)).toBe('ready');
    expect(loader.getState(2, 0)).toBe('none');
    expect(list.targetMapping!.get(1)).toBe(loader.getTarget(1));
    expect(view.cell!.obj!.description).toBe(`2 of 8 tiles · ${list.count} particles`);
  });

  it('keeps a tile-ordered pack opened from a file in a session, read from it by range when opened again', async () => {
    const bytes = writeParticlePack(createList(), { tiles: 2 });
    const file = Asset.File(new File([bytes], 'lattice.mspk'));
    const tiles = await plugin
      .build()
      .toRoot()
      .apply(ParticleTileSetFromPack, { source: { name: 'file', params: file } } as never)
      .commit();
    expect(tiles.data!.set.count).toBe(Count);
    // the file held (not read): the session's zip keeps it
    expect(plugin.managers.asset.get(file)?.file.size).toBe(bytes.byteLength);
    const session = await plugin.managers.snapshot.serialize({
      type: 'molx',
      params: { camera: false, canvas3d: false },
    });

    const other = new PluginContext({ behaviors: [], registry: [...DefaultRegistry, ...ParticlesRegistry] });
    await other.init();
    try {
      await other.managers.snapshot.open(new File([session], 'session.molx'));
      const cells = other.state.data.select(StateSelection.Generators.ofTransformer(ParticleTileSetFromPack));
      expect(cells.length).toBe(1);
      expect(cells[0].status).toBe('ok');
      expect((cells[0].obj!.data as typeof tiles.data)!.set.count).toBe(Count);
    } finally {
      other.dispose();
    }
    // the set removed: its file no longer held
    await plugin.build().delete(tiles.ref).commit();
    expect(plugin.managers.asset.has(file)).toBe(false);
  });

  it("picks the grid by default: a tile-ordered list's written one, else about a page of particles a tile", async () => {
    const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
    const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, {}).commit();
    const expected = ParticleTiles.getAutoCells(list.data!);
    expect(tiles.data!.set.grid.dimensions).toEqual([expected, expected, expected]);
    expect(tiles.cell!.obj!.description).toContain(`${expected}³ grid`);
  });

  it('keeps the type and theme of a representation of a slot view of no tiles yet, drawn as its tiles come', async () => {
    const hierarchical: ParticleList = {
      ...createList(),
      entities: new Int32Array(Count),
      entityInfo: new Map([[0, { name: 'entity' }]]),
      compartments: new Int32Array(Count),
      compartmentInfo: new Map([[0, { name: 'root.cyto' }]]),
    };
    const list = await plugin.build().toRoot().apply(ListOf, { list: hierarchical }).commit();
    const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();
    // none of the tiles yet (as the camera seeing none of them): a slot list of no pages
    const view = await plugin
      .build()
      .to(tiles.ref)
      .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids: '' } }, frame: 0, slots: true })
      .commit();
    expect(view.data!.pages!.pages.length).toBe(0);
    const repr = await plugin
      .build()
      .to(view.ref)
      .apply(ParticlesRepresentation3D, {
        type: { name: 'target', params: {} },
        colorTheme: { name: 'particle-hierarchy', params: {} },
      } as never)
      .commit();
    const params = () => repr.cell!.transform.params as { type: { name: string }; colorTheme: { name: string } };
    expect(params().type.name).toBe('target');
    expect(params().colorTheme.name).toBe('particle-hierarchy');

    await plugin
      .build()
      .to(view.ref)
      .update(ParticleTilesView, (old) => ({ ...old, tiles: { name: 'all', params: {} } }))
      .commit();
    expect(params().type.name).toBe('target');
    expect(params().colorTheme.name).toBe('particle-hierarchy');
    expect(repr.data!.repr.renderObjects.length).toBeGreaterThan(0);
  });

  it('updates the view as its tiles change, the targets drawn once whichever tiles', async () => {
    const tiles = await tileSet({ delay: 5 });
    const view = await plugin
      .build()
      .to(tiles.ref)
      .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids: '0' } }, frame: 0, slots: false })
      .commit();
    const repr = await plugin
      .build()
      .to(view.ref)
      .apply(ParticlesRepresentation3D, { type: { name: 'target', params: {} } })
      .commit();
    const visuals = () => repr.data!.repr.renderObjects.length;
    const one = visuals();
    expect(one).toBeGreaterThan(0);

    await plugin
      .build()
      .to(view.ref)
      .update(ParticleTilesView, (old) => ({ ...old, tiles: { name: 'all', params: {} } }))
      .commit();
    expect(view.cell!.obj!.data.count).toBe(Count);
    // a visual a target, the same count for the whole as for a tile of both targets
    expect(visuals()).toBe(one);
    // by a box: the tiles of x < 45
    await plugin
      .build()
      .to(view.ref)
      .update(ParticleTilesView, (old) => ({
        ...old,
        tiles: { name: 'box', params: { min: Vec3.create(0, 0, 0), max: Vec3.create(40, 100, 100) } },
      }))
      .commit();
    const v = ParticleTileView.of(view.cell!.obj!.data)!;
    expect(v.tiles.length).toBe(4);
    const p = Vec3();
    for (let i = 0; i < v.list.count; ++i) expect(GridPositions.get(v.list.coordinates, i, p)[0]).toBeLessThan(45);
  });

  it("keeps a target's instances where a view of other tiles has its part the same, their indices the new ones", async () => {
    // targets by layer: z below 45 (tiles 0-3) target 1, above (tiles 4-7) target 2
    const base = createList();
    const targets = Int32Array.from({ length: Count }, (_, i) => (Math.floor(i / 16) < 2 ? 1 : 2));
    const l = await plugin
      .build()
      .toRoot()
      .apply(ListOf, { list: { ...base, targets } })
      .commit();
    const tiles = await plugin.build().to(l.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();
    const view = await plugin
      .build()
      .to(tiles.ref)
      .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids: '4' } }, frame: 0, slots: false })
      .commit();
    const repr = await plugin
      .build()
      .to(view.ref)
      .apply(ParticlesRepresentation3D, { type: { name: 'target', params: { renderMerged: false } } } as any)
      .commit();
    const r = repr.data!.repr;
    const ro = r.renderObjects[0];
    const version = () => (ro.values as any).aTransform.ref.version;
    const before = version();
    const setIds = (ids: string) =>
      plugin
        .build()
        .to(view.ref)
        .update(ParticleTilesView, (old) => ({ ...old, tiles: { name: 'ids', params: { ids } } }))
        .commit();

    // tile 0 of target 1 added: target 2's part the same (tile 4), after target 1's particles in the list now
    await setIds('0, 4');
    expect(r.renderObjects).toHaveLength(2);
    expect(r.renderObjects).toContain(ro);
    expect(version()).toBe(before);
    const list = view.cell!.obj!.data;
    const loci = r.getLoci({ objectId: ro.id, instanceId: 0, groupId: 0 }) as Particle.Loci;
    expect(Particle.isLoci(loci) && loci.particles).toBe(list);
    const first = OrderedSet.getAt(loci.indices, 0);
    expect(ParticleValues.get(list.targets, first)).toBe(2);
    expect(first).toBeGreaterThanOrEqual(tiles.data!.set.tiles.get(0)!.count);

    // tile 5 of target 2 added: its part changed, its instances built again
    await setIds('0, 4, 5');
    expect(version()).not.toBe(before);
  });

  it('colors a particle by its index in the whole in views of tiles, so the same in views of other tiles', async () => {
    const tiles = await tileSet();
    const { set, loader } = tiles.data!;
    loader.setWanted([0, 1].map((id) => ({ id, frame: 0 })));
    await loader.whenIdle();
    const both = ParticleTileView.create(set, loader, [0, 1], 0).list;
    const one = ParticleTileView.create(set, loader, [1], 0).list;
    const colorOf = (list: ParticleList) => {
      const theme = ParticleIndexColorThemeProvider.factory(
        { particles: list },
        PD.getDefaultValues(ParticleIndexColorThemeProvider.getParams({ particles: list })),
      );
      return (i: number) => theme.color(Particle.Location(list, i), false);
    };
    const inBoth = new Map<number, number>();
    const c2 = colorOf(both);
    for (let i = 0; i < both.count; ++i) inBoth.set(both.keys![i], c2(i));
    const c1 = colorOf(one);
    for (let i = 0; i < one.count; ++i) expect(c1(i)).toBe(inBoth.get(one.keys![i]));
    // a list not a view: by its index
    const list = createList();
    const c = colorOf(list);
    expect(c(0)).not.toBe(c(1));
  });

  it('steps a trajectory: the same particles posed otherwise', async () => {
    const trajectory = new ArrayParticleTrajectory([createList(), createList(7)]);
    const t = await plugin.build().toRoot().apply(TrajectoryOf, { trajectory }).commit();
    const tiles = await plugin.build().to(t.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();
    expect(tiles.data!.set.frameCount).toBe(2);
    const view = await plugin
      .build()
      .to(tiles.ref)
      .apply(ParticleTilesView, { tiles: { name: 'all', params: {} }, frame: 0, slots: false })
      .commit();
    const first = view.data!;
    await plugin
      .build()
      .to(view.ref)
      .update(ParticleTilesView, (old) => ({ ...old, frame: 1 }))
      .commit();
    const second = view.cell!.obj!.data;
    expect(second).not.toBe(first);
    expect(Particle.isSamePoses(first, second)).toBe(true);
    expect(view.cell!.obj!.description).toContain('frame 2 of 2');
    const p = Vec3();
    const q = Vec3();
    expect(GridPositions.get(second.coordinates, 0, q)[0] - GridPositions.get(first.coordinates, 0, p)[0]).toBeCloseTo(
      7,
      4,
    );
  });

  it('caps the particles a camera view picks: the tiles by priority while within the cap, at least one', async () => {
    const tiles = await tileSet({ cells: 4 });
    const set = tiles.data!.set;
    const camera = {
      snapshot: {
        mode: 'perspective' as const,
        fov: Math.PI / 4,
        position: Vec3.create(45, 45, 400),
        target: Vec3.create(45, 45, 0),
        up: Vec3.create(0, 1, 0),
        minNear: 1,
      },
      width: 512,
      height: 512,
    };
    const pick = (maxParticles: number) =>
      pickParticleTiles(
        set,
        {
          name: 'camera',
          params: {
            minPixels: 0,
            occlusion: false,
            opticalDepth: 7,
            follow: false,
            budget: 0,
            maxParticles,
            stepMs: 0,
            camera,
          },
        },
        plugin,
      );
    const all = pick(0);
    expect(all.length).toBeGreaterThan(4);
    const countOf = (ids: number[]) => ids.reduce((n, id) => n + set.tiles.get(id)!.count, 0);
    const cap = Math.floor(countOf(all) / 2);
    const capped = pick(cap);
    // the uncapped pick's first tiles, as many as fit
    expect(capped).toEqual(all.slice(0, capped.length));
    expect(countOf(capped)).toBeLessThanOrEqual(cap);
    expect(countOf(all.slice(0, capped.length + 1))).toBeGreaterThan(cap);
    // a cap below a tile: the first alone
    expect(pick(1)).toEqual(all.slice(0, 1));
  });

  it("picks again where the camera's near plane alone moved: a pick by the camera's minNear too", async () => {
    const tiles = await tileSet();
    const set = tiles.data!.set;
    const pick = (minNear: number) =>
      pickParticleTiles(
        set,
        {
          name: 'camera',
          params: {
            minPixels: 0,
            occlusion: false,
            opticalDepth: 7,
            follow: false,
            budget: 0,
            maxParticles: 0,
            stepMs: 0,
            camera: {
              snapshot: {
                mode: 'perspective' as const,
                fov: Math.PI / 4,
                position: Vec3.create(45, 45, 400),
                target: Vec3.create(45, 45, 0),
                up: Vec3.create(0, 1, 0),
                minNear,
              },
              width: 512,
              height: 512,
            },
          },
        },
        plugin,
      );
    const near = pick(1);
    expect(near.length).toBeGreaterThan(0);
    // the near plane past the whole lattice (z 0-90, the camera at 400): none of its tiles in front of it
    const far = pick(1000);
    expect(far).not.toEqual(near);
    // views of the set by these two cameras in turn: each its own pick kept, not made again
    expect(pick(1)).toBe(near);
    expect(pick(1000)).toBe(far);
    expect(pick(1)).toBe(near);
  });

  it("picks a set's coarse levels: a cut of the octree, nodes where small on screen, coarser within the cap", async () => {
    // the lattice a tile a particle (4³ tiles), its levels: 8 nodes of 8 tiles, the root of all
    const bytes = writeParticlePack(createList(), { tiles: 4, levels: true });
    const set = await createParticleTileSetFromPack(await openParticlePackStream(ParticlePackSource.fromBytes(bytes)));
    const levels = set.levels!;
    const camera = {
      snapshot: {
        mode: 'perspective' as const,
        fov: Math.PI / 4,
        position: Vec3.create(45, 45, 400),
        target: Vec3.create(45, 45, 0),
        up: Vec3.create(0, 1, 0),
        minNear: 1,
      },
      width: 512,
      height: 512,
    };
    const pick = (detailPixels: number, maxParticles = 0) =>
      pickParticleTiles(
        set,
        {
          name: 'camera',
          params: {
            minPixels: 0,
            occlusion: false,
            opticalDepth: 7,
            follow: false,
            budget: 0,
            detailPixels,
            maxParticles,
            stepMs: 0,
            camera,
          },
        },
        plugin,
      );
    const countOf = (ids: number[]) => ids.reduce((n, id) => n + set.tiles.get(id)!.count, 0);
    // the default detail, the tiles large enough on screen: the visible tiles themselves
    const visible = pick(4);
    expect(visible.length).toBe(64);
    expect(visible.every((id) => levels.base.tiles.has(id))).toBe(true);
    // each visible tile shown once: by itself or by one node over it
    const showsEachOnce = (ids: number[]) => {
      const picked = new Set(ids);
      for (const tile of visible) {
        let n = picked.has(tile) ? 1 : 0;
        for (let level = 1; level <= levels.count; ++level) {
          const node = levels.getNode(level, ParticleTiles.getLevelCell(4, tile, level));
          if (node !== undefined && picked.has(node)) n += 1;
        }
        expect(n, `tile ${tile}`).toBe(1);
      }
    };
    // small on screen: the root alone
    expect(pick(1e9)).toEqual([levels.getNode(2, 0)]);
    // within a cap: coarser (a node of 8 tiles a particle) until it fits, each tile still shown
    const capped = pick(4, 8);
    expect(countOf(capped)).toBeLessThanOrEqual(8);
    expect(capped.every((id) => !levels.base.tiles.has(id))).toBe(true);
    showsEachOnce(capped);
    // a cap below that coarsest the detail lowered at most: the least important nodes left out, at least one
    const less = pick(4, 3);
    expect(countOf(less)).toBeLessThanOrEqual(3);
    expect(less.length).toBeGreaterThan(0);
  });

  it('keeps the tiles another view of the set wants or shows: a capped camera view lets go of none of them', async () => {
    const tiles = await tileSet();
    const { set, loader } = tiles.data!;
    // a camera view on a corner, capped below what is loaded: it lets go of the tiles it neither wants nor shows
    const camera = (x: number) => ({
      snapshot: {
        mode: 'perspective' as const,
        fov: Math.PI / 64,
        position: Vec3.create(x, 15, 400),
        target: Vec3.create(x, 15, 0),
        up: Vec3.create(0, 1, 0),
        minNear: 1,
      },
      width: 512,
      height: 512,
    });
    const capped = (x: number) => ({
      tiles: {
        name: 'camera' as const,
        params: { minPixels: 0, occlusion: false, maxParticles: 1, camera: camera(x) },
      },
      frame: 0,
    });
    const byCamera = await plugin.build().to(tiles.ref).apply(ParticleTilesView, capped(15)).commit();
    const seen = new Set(ParticleTileView.of(byCamera.data!)!.tiles);
    // a view by id of the tiles the camera view does not show
    const others = [...set.tiles.keys()].filter((id) => !seen.has(id)).sort((a, b) => a - b);
    expect(others.length).toBeGreaterThan(0);
    const byId = await plugin
      .build()
      .to(tiles.ref)
      .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids: others.join(', ') } }, frame: 0 })
      .commit();
    // the camera view's next step: the other view's tiles kept
    await plugin
      .build()
      .to(byCamera.ref)
      .update(ParticleTilesView, (old) => ({ ...old, ...capped(16) }))
      .commit();
    for (const id of others) expect(loader.getState(id, 0)).toBe('ready');
    expect([...ParticleTileView.of(byId.cell!.obj!.data)!.tiles]).toEqual(others);
    // the view by id gone: its tiles let go at the camera view's next step
    await plugin.build().delete(byId.ref).commit();
    await plugin
      .build()
      .to(byCamera.ref)
      .update(ParticleTilesView, (old) => ({ ...old, ...capped(15) }))
      .commit();
    for (const id of others) expect(loader.getState(id, 0)).toBe('none');
  });

  it('picks the tiles the camera sees, others as it moves', async () => {
    const tiles = await tileSet();
    // a narrow view along -z: the lattice's corner at x, y 0-45, then its other one at 45-90
    const camera = (x: number, y: number) => ({
      snapshot: {
        mode: 'perspective' as const,
        fov: Math.PI / 64,
        position: Vec3.create(x, y, 400),
        target: Vec3.create(x, y, 0),
        up: Vec3.create(0, 1, 0),
        minNear: 1,
      },
      width: 512,
      height: 512,
    });
    const view = await plugin
      .build()
      .to(tiles.ref)
      .apply(ParticleTilesView, {
        tiles: { name: 'camera', params: { minPixels: 0, occlusion: true, camera: camera(15, 15) } },
        frame: 0,
      })
      .commit();
    const ids = () => [...ParticleTileView.of(view.cell!.obj!.data)!.tiles];
    const set = tiles.data!.set;
    const xyOf = (id: number) => {
      const t = set.tiles.get(id)!;
      return [t.box.min[0] < 45 ? 0 : 1, t.box.min[1] < 45 ? 0 : 1];
    };
    expect(
      ids()
        .map(xyOf)
        .every(([x, y]) => x === 0 && y === 0),
    ).toBe(true);
    await plugin
      .build()
      .to(view.ref)
      .update(ParticleTilesView, (old) => ({
        ...old,
        tiles: { name: 'camera', params: { minPixels: 0, occlusion: true, camera: camera(75, 75) } },
      }))
      .commit();
    expect(ids().length).toBeGreaterThan(0);
    expect(
      ids()
        .map(xyOf)
        .every(([x, y]) => x === 1 && y === 1),
    ).toBe(true);
  });

  it('hides the tiles behind dense ones where it computed their opacities', async () => {
    // a dense lattice: 8000 carbons a particle (a sphere of 31 Å) 30 Å apart
    const base = createList();
    const dense = (n: number): ParticleTarget => ({
      kind: 'sites',
      sites: {
        count: n,
        positions: new Float32Array(n * 3),
        radii: new Float32Array(n).fill(2),
        elements: new Uint8Array(n).fill(6),
      },
    });
    const list = {
      ...base,
      targetMapping: new Map([
        [1, dense(8000)],
        [2, dense(8000)],
      ]),
    };
    const l = await plugin.build().toRoot().apply(ListOf, { list }).commit();
    const tiles = await plugin.build().to(l.ref).apply(ParticleTileSetFromList, { cells: 4, opacity: true }).commit();
    expect([...tiles.data!.set.tiles.values()].every((t) => t.opacity)).toBe(true);
    const view = (occlusion: boolean) => ({
      tiles: {
        name: 'camera' as const,
        params: {
          minPixels: 0,
          occlusion,
          camera: {
            snapshot: {
              mode: 'perspective' as const,
              // narrow: every ray through the front tile (the one behind's reach, 2 Å past its box, else seen past it)
              fov: Math.PI / 512,
              // along a column of molecules (the tiles' corner at 45, 45 between them)
              position: Vec3.create(30, 30, 1000),
              target: Vec3.create(30, 30, 0),
              up: Vec3.create(0, 1, 0),
              minNear: 1,
            },
            width: 512,
            height: 512,
          },
        },
      },
      frame: 0,
    });
    const without = await plugin.build().to(tiles.ref).apply(ParticleTilesView, view(false)).commit();
    const withOcclusion = await plugin.build().to(tiles.ref).apply(ParticleTilesView, view(true)).commit();
    const n = (c: typeof without) => ParticleTileView.of(c.data!)!.tiles.length;
    // a column of four tiles of 22.5 Å in view, each a molecule's 0.27 an Å: the front two deep enough together,
    // the third's reach (2 Å past its box) within the second's cell, so kept; the last hidden
    expect(n(without)).toBe(4);
    expect(n(withOcclusion)).toBe(3);
  });
});

describe('particle tile views following the camera', () => {
  // a canvas of a camera alone: what the behavior and the camera mode read of it
  function stubCanvas(plugin: PluginContext) {
    let snapshot = {
      mode: 'perspective' as const,
      fov: Math.PI / 64,
      position: Vec3.create(15, 15, 400),
      target: Vec3.create(15, 15, 0),
      up: Vec3.create(0, 1, 0),
      minNear: 1,
    };
    const changed = new Subject<void>();
    const camera = { changed, getSnapshot: () => snapshot, viewport: { width: 512, height: 512 } };
    (plugin as unknown as { canvas3d: unknown }).canvas3d = { camera, dispose: () => {} };
    (plugin as unknown as { canvas3dInit: Subject<boolean> }).canvas3dInit.next(true);
    return {
      lookAt(x: number, y: number) {
        snapshot = { ...snapshot, position: Vec3.create(x, y, 400), target: Vec3.create(x, y, 0) };
        changed.next();
      },
    };
  }
  const settled = async (done: () => boolean) => {
    for (let i = 0; i < 100 && !done(); ++i) await new Promise((r) => setTimeout(r, 10));
  };

  it('picks the tiles of a following view again as the camera moves, a view of a camera given as it is', async () => {
    const plugin = new PluginContext({
      behaviors: [PluginSpec.Behavior(ParticleTilesFollowCamera, { intervalMs: 0 })],
      registry: [...DefaultRegistry, ...ParticlesRegistry],
    });
    await plugin.init();
    try {
      const canvas = stubCanvas(plugin);
      const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();
      const set = tiles.data!.set;
      const corner = (id: number) => {
        const t = set.tiles.get(id)!;
        return `${t.box.min[0] < 45 ? 0 : 1}${t.box.min[1] < 45 ? 0 : 1}`;
      };
      const view = (follow: boolean, given: boolean) =>
        plugin
          .build()
          .to(tiles.ref)
          .apply(ParticleTilesView, {
            tiles: {
              name: 'camera',
              params: {
                minPixels: 0,
                occlusion: true,
                follow,
                camera: given
                  ? { snapshot: plugin.canvas3d!.camera.getSnapshot(), width: 512, height: 512 }
                  : undefined,
              },
            },
            frame: 0,
          })
          .commit();
      const following = await view(true, false);
      const fixed = await view(true, true);
      const unfollowed = await view(false, false);
      const corners = (v: typeof following) => [
        ...new Set([...ParticleTileView.of(v.cell!.obj!.data)!.tiles].map(corner)),
      ];
      expect(corners(following)).toEqual(['00']);

      // a pick that keeps the tiles shown and brings others: those loaded (the loader given the new pick), the view
      // stepping to them
      canvas.lookAt(45, 45);
      await settled(() => corners(following).length === 4);
      expect(corners(following).sort()).toEqual(['00', '01', '10', '11']);

      canvas.lookAt(75, 75);
      await settled(() => corners(following).join() === '11');
      expect(corners(following)).toEqual(['11']);
      // a camera given, or following off: as made
      expect(corners(fixed)).toEqual(['00']);
      expect(corners(unfollowed)).toEqual(['00']);
      // at rest where it was: no update
      const version = following.cell!.transform.version;
      canvas.lookAt(75, 75);
      await new Promise((r) => setTimeout(r, 50));
      expect(following.cell!.transform.version).toBe(version);
    } finally {
      plugin.dispose();
    }
  });

  it("lets go of a tile set's loader as the set goes: none of its subscriptions kept", async () => {
    const plugin = new PluginContext({
      behaviors: [PluginSpec.Behavior(ParticleTilesFollowCamera, { intervalMs: 0 })],
      registry: [...DefaultRegistry, ...ParticlesRegistry],
    });
    // the loader's live subscriptions counted, from its set's creation (before the behavior's, which follows it as
    // it is made: subscribed before init, which registers the behavior)
    let live = 0;
    plugin.state.data.events.object.created.subscribe(({ obj }) => {
      if (!SO.Particle.TileSet.is(obj)) return;
      const loader = obj.data.loader;
      const subscribe = loader.subscribe;
      loader.subscribe = (listener) => {
        const unsubscribe = subscribe(listener);
        live += 1;
        return () => {
          live -= 1;
          return unsubscribe();
        };
      };
    });
    await plugin.init();
    try {
      const canvas = stubCanvas(plugin);
      const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();
      await plugin

        .build()
        .to(tiles.ref)
        .apply(ParticleTilesView, {
          tiles: { name: 'camera', params: { minPixels: 0, occlusion: false, follow: true, camera: undefined } },
          frame: 0,
        })
        .commit();
      canvas.lookAt(45, 45);
      await settled(() => live > 0);
      expect(live).toBe(1);
      await plugin.build().delete(tiles.ref).commit();
      expect(live).toBe(0);
    } finally {
      plugin.dispose();
    }
  });

  it('adds the tiles a following view picks in steps of its budget, as they load', async () => {
    const plugin = new PluginContext({
      behaviors: [PluginSpec.Behavior(ParticleTilesFollowCamera, { intervalMs: 0 })],
      registry: [...DefaultRegistry, ...ParticlesRegistry],
    });
    await plugin.init();
    try {
      const canvas = stubCanvas(plugin);
      // a wide view: all tiles in it
      canvas.lookAt(45, 45);
      const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();
      const params = {
        tiles: {
          name: 'camera' as const,
          params: { minPixels: 0, occlusion: false, follow: true, budget: 1, camera: undefined },
        },
        frame: 0,
      };
      const picked = pickParticleTiles(tiles.data!.set, params.tiles, plugin);
      expect(picked.length).toBeGreaterThan(2);
      const view = await plugin.build().to(tiles.ref).apply(ParticleTilesView, params).commit();
      const shown = () => ParticleTileView.of(view.cell!.obj!.data)!.tiles.length;
      // a tile a step (the budget below a tile's particles: one), the first at once
      expect(shown()).toBe(1);
      // each step's count as the view's object is updated (not polled: a fast step between two polls is missed)
      const seen = new Set<number>([shown()]);
      const steps = plugin.state.data.events.object.updated.subscribe((e) => {
        if (e.ref === view.ref) seen.add(ParticleTileView.of(e.obj.data)!.tiles.length);
      });
      try {
        await settled(() => shown() === picked.length);
      } finally {
        steps.unsubscribe();
      }
      expect(shown()).toBe(picked.length);
      // a step at a time: the counts between seen
      expect(seen.size).toBeGreaterThan(2);
    } finally {
      plugin.dispose();
    }
  });

  it('steps a following view of a set made again (a new loader) with the camera still', async () => {
    const plugin = new PluginContext({
      behaviors: [PluginSpec.Behavior(ParticleTilesFollowCamera, { intervalMs: 0 })],
      registry: [...DefaultRegistry, ...ParticlesRegistry],
    });
    await plugin.init();
    try {
      const canvas = stubCanvas(plugin);
      canvas.lookAt(45, 45);
      const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2, delay: 10 }).commit();
      const params = {
        tiles: {
          name: 'camera' as const,
          params: { minPixels: 0, occlusion: false, follow: true, budget: 1, camera: undefined },
        },
        frame: 0,
      };
      const picked = pickParticleTiles(tiles.data!.set, params.tiles, plugin).length;
      expect(picked).toBeGreaterThan(2);
      const view = await plugin.build().to(tiles.ref).apply(ParticleTilesView, params).commit();
      const shown = () => ParticleTileView.of(view.cell!.obj!.data)!.tiles.length;
      await settled(() => shown() === picked);
      expect(shown()).toBe(picked);
      // the set of other params (it has no update: made again, a new loader), the view its first step of it
      await plugin
        .build()
        .to(tiles.ref)
        .update(ParticleTileSetFromList, (p) => ({ ...p, delay: 20 }))
        .commit();
      expect(shown()).toBeLessThan(picked);
      await settled(() => shown() === picked);
      expect(shown()).toBe(picked);
    } finally {
      plugin.dispose();
    }
  });

  it('keeps a following view as it is where its next pick is all new and none of it ready: never empty', async () => {
    const plugin = new PluginContext({
      behaviors: [PluginSpec.Behavior(ParticleTilesFollowCamera, { intervalMs: 0 })],
      registry: [...DefaultRegistry, ...ParticlesRegistry],
    });
    await plugin.init();
    try {
      const canvas = stubCanvas(plugin);
      const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
      // tiles served after a while, as a fetch would
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2, delay: 40 }).commit();
      const view = await plugin
        .build()
        .to(tiles.ref)
        .apply(ParticleTilesView, {
          tiles: { name: 'camera', params: { minPixels: 0, occlusion: true, follow: true, camera: undefined } },
          frame: 0,
        })
        .commit();
      const shown = () => [...ParticleTileView.of(view.cell!.obj!.data)!.tiles];
      const first = shown();
      expect(first.length).toBeGreaterThan(0);
      // the view's updates counted: none while it waits (the same tiles again every frame), its steps after
      const key = () => shown().join(',');
      let same = 0;
      let last = key();
      const updates = plugin.state.data.events.object.updated.subscribe((e) => {
        if (e.ref !== view.ref) return;
        if (key() === last) same += 1;
        last = key();
      });
      // a pick of other tiles alone: the view the first's until they are ready, then theirs
      canvas.lookAt(75, 75);
      let least = first.length;
      try {
        await settled(() => {
          least = Math.min(least, shown().length);
          return shown().some((id) => !first.includes(id));
        });
      } finally {
        updates.unsubscribe();
      }
      expect(least).toBeGreaterThan(0);
      expect(shown().some((id) => !first.includes(id))).toBe(true);
      expect(same).toBe(0);
    } finally {
      plugin.dispose();
    }
  });

  it('lets go of the tiles a capped following view passed, past twice its cap', async () => {
    const plugin = new PluginContext({
      behaviors: [PluginSpec.Behavior(ParticleTilesFollowCamera, { intervalMs: 0 })],
      registry: [...DefaultRegistry, ...ParticlesRegistry],
    });
    await plugin.init();
    try {
      const canvas = stubCanvas(plugin);
      const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 4 }).commit();
      const loader = tiles.data!.loader;
      const cap = 2;
      const view = await plugin
        .build()
        .to(tiles.ref)
        .apply(ParticleTilesView, {
          tiles: {
            name: 'camera',
            params: { minPixels: 0, occlusion: false, follow: true, maxParticles: cap, camera: undefined },
          },
          frame: 0,
        })
        .commit();
      const shown = () => [...ParticleTileView.of(view.cell!.obj!.data)!.tiles];
      const passed = new Set<number>(shown());
      for (const [x, y] of [
        [75, 15],
        [75, 75],
        [15, 75],
        [45, 45],
      ]) {
        const before = shown().join();
        canvas.lookAt(x, y);
        await settled(() => shown().join() !== before);
        for (const id of shown()) passed.add(id);
      }
      // more passed than kept: the loader within twice the cap and what the view shows and wants
      expect(passed.size).toBeGreaterThan(2 * cap);
      expect(loader.getLoaded().particles).toBeLessThanOrEqual(2 * cap + shown().length + cap);
    } finally {
      plugin.dispose();
    }
  });

  it("steps a following view by time: a step's particles its set's (by the last steps' time), with stepMs, at most its budget", async () => {
    // the particles of a step of about stepMs at a rate (particles a ms), at least one; none without a rate or a time
    expect(getParticleTilesStepParticles(2000, 10)).toBe(20_000);
    expect(getParticleTilesStepParticles(0.01, 10)).toBe(1);
    expect(getParticleTilesStepParticles(undefined, 10)).toBeUndefined();
    expect(getParticleTilesStepParticles(2000, 0)).toBeUndefined();

    const plugin = new PluginContext({ behaviors: [], registry: [...DefaultRegistry, ...ParticlesRegistry] });
    await plugin.init();
    try {
      const canvas = stubCanvas(plugin);
      canvas.lookAt(45, 45);
      const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();
      const set = tiles.data!.set;
      const params = (stepMs: number) => ({
        tiles: {
          name: 'camera' as const,
          params: { minPixels: 0, occlusion: false, follow: true, budget: 1_000_000, stepMs, camera: undefined },
        },
        frame: 0,
      });
      const picked = pickParticleTiles(set, params(10).tiles, plugin);
      expect(picked.length).toBeGreaterThan(2);
      const shown = (v: { cell?: { obj?: { data: unknown } } }) =>
        ParticleTileView.of(v.cell!.obj!.data as never)!.tiles.length;
      // a step of a particle (the last steps' time): a tile a step, the budget not reached
      setParticleTilesStep(set, 1);
      const byTime = await plugin.build().to(tiles.ref).apply(ParticleTilesView, params(10)).commit();
      expect(shown(byTime)).toBe(1);
      // stepMs 0: the budget alone, all ready at once
      const byBudget = await plugin.build().to(tiles.ref).apply(ParticleTilesView, params(0)).commit();
      expect(shown(byBudget)).toBe(picked.length);
      // none set: the budget
      setParticleTilesStep(set, undefined);
      const unset = await plugin.build().to(tiles.ref).apply(ParticleTilesView, params(10)).commit();
      expect(shown(unset)).toBe(picked.length);
    } finally {
      plugin.dispose();
    }
  });

  it('leaves out the tiles behind others by the occlusion threshold (opticalDepth): fewer at a lower one', async () => {
    const plugin = new PluginContext({ behaviors: [], registry: [...DefaultRegistry, ...ParticlesRegistry] });
    await plugin.init();
    try {
      const canvas = stubCanvas(plugin);
      canvas.lookAt(45, 45);
      // four layers of tiles along the view (z), each tile an optical depth of 2 through it
      const base = ParticleTiles.fromList(createList(), [2, 2, 4]);
      const depth = base.grid.cellSize[2];
      const tiles = new Map(
        [...base.tiles].map(([id, t]) => [
          id,
          { ...t, opacity: { dimensions: [1, 1, 1] as const, density: Float32Array.of(2 / depth) } },
        ]),
      );
      const set = { ...base, tiles };
      const pick = (opticalDepth: number) =>
        pickParticleTiles(
          set,
          {
            name: 'camera',
            params: {
              minPixels: 0,
              occlusion: true,
              opticalDepth,
              follow: false,
              budget: 0,
              stepMs: 0,
              camera: undefined,
            },
          },
          plugin,
        );
      // a ray through all four 8: not stopped at 7; at 1 by the front layer, the two behind the next left out (the next
      // within the reach of its targets)
      expect(pick(7).length).toBe(16);
      expect(pick(1).length).toBe(8);
      expect(pick(1).every((id) => set.tiles.get(id)!.box.max[2] > 50)).toBe(true);
    } finally {
      plugin.dispose();
    }
  });
});
