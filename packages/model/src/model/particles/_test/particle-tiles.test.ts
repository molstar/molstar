/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { Box3D } from '@molstar/core/math/geometry/primitives/box3d';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Task } from '@molstar/core/task';
import { CustomProperties } from '../../custom-property.js';
import {
  Particle,
  type ParticleList,
  ParticleRotations,
  ParticleSites,
  type ParticleTarget,
  ParticleValues,
} from '../particle-list.js';
import { ParticleTileCatalog, ParticleTileGrid, ParticleTiles } from '../particle-tiles.js';
import { ArrayParticleTrajectory } from '../particle-trajectory.js';

const Count = 200;

/** a deterministic spread of `Count` particles over 0-100 Å, targets 1 and 2 in runs of ten, an attribute */
function createList(shift = 0): ParticleList {
  const coordinates = new Float32Array(Count * 3);
  const rotations = new Float32Array(Count * 4);
  for (let i = 0; i < Count; ++i) {
    coordinates.set([((i * 37) % 101) + shift, (i * 53) % 97, (i * 11) % 89], i * 3);
    rotations.set([0, 0, Math.sin(i), Math.cos(i)], i * 4);
  }
  const targets = Int32Array.from({ length: Count }, (_, i) => (Math.floor(i / 10) % 2 ? 2 : 1));
  const sites = (n: number): ParticleTarget => ({
    kind: 'sites',
    sites: { count: n, positions: new Float32Array(n * 3), radii: new Float32Array(n).fill(1.5) },
  });
  return {
    label: 'test',
    count: Count,
    targets,
    targetInfo: new Map([
      [1, { name: 'one', entity: 0 }],
      [2, { name: 'two', entity: 1 }],
    ]),
    targetMapping: new Map([
      [1, sites(3)],
      [2, sites(5)],
    ]),
    entities: ParticleValues.ofRuns(Int32Array.of(0, 100, Count), Int32Array.of(0, 1)),
    entityInfo: new Map([
      [0, { name: 'e0' }],
      [1, { name: 'e1' }],
    ]),
    coordinates,
    rotations,
    attributes: new Map([
      [
        'score',
        {
          column: Column.ofFloatArray(Float32Array.from({ length: Count }, (_, i) => i / 10)),
          label: 'Score',
          min: 0,
          max: 20,
        },
      ],
    ]),
    getParticleLabel: (index) => `particle ${index}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

/** the index in `list` of a tile's particle `i`, by its label (`Particle.subset` labels a particle by its index) */
const sourceIndex = (tile: ParticleList, i: number) => Number(tile.getParticleLabel(i).split(' ')[1]);

describe('Particle.subset', () => {
  it('gathers the particles of indices in their order, their values, the targets the same objects', () => {
    const list = createList();
    const indices = Int32Array.of(150, 3, 4, 120);
    const sub = Particle.subset(list, indices, { label: 'part', keys: Int32Array.of(7, 8, 9, 10) });
    expect(sub.count).toBe(4);
    expect(sub.label).toBe('part');
    expect(Array.from(sub.keys!)).toEqual([7, 8, 9, 10]);
    expect(sub.targetInfo).toBe(list.targetInfo);
    expect(sub.targetMapping).toBe(list.targetMapping);
    expect(sub.entityInfo).toBe(list.entityInfo);
    const p = Vec3();
    const q = Vec3();
    indices.forEach((src, i) => {
      expect(Array.from(GridPositions.get(sub.coordinates, i, p))).toEqual(
        Array.from(GridPositions.get(list.coordinates, src, q)),
      );
      expect(Array.from(sub.rotations!.subarray(i * 4, i * 4 + 4))).toEqual(
        Array.from(list.rotations!.subarray(src * 4, src * 4 + 4)),
      );
      expect(ParticleValues.get(sub.targets, i)).toBe(ParticleValues.get(list.targets, src));
      expect(ParticleValues.get(sub.entities!, i)).toBe(ParticleValues.get(list.entities!, src));
      expect(sub.attributes!.get('score')!.column.value(i)).toBeCloseTo(src / 10, 5);
      expect(sub.getParticleLabel(i)).toBe(`particle ${src}`);
    });
    // runs where fewer than half: the first two of entity 1, the next two of 0
    expect(ParticleValues.isRuns(sub.entities!)).toBe(false);
    const runs = Particle.subset(list, Int32Array.of(0, 1, 2, 3, 4, 5));
    expect(ParticleValues.isRuns(runs.targets)).toBe(true);
    expect(sub.fibers).toBeUndefined();
  });

  it('keeps positions on a grid exact and rotations as snorm16', () => {
    const base = createList();
    const list = {
      ...base,
      coordinates: GridPositions.fromFloat(base.coordinates as Float32Array),
      rotations: ParticleRotations.toSnorm16(base.rotations!),
    };
    const indices = Int32Array.of(9, 0, 199);
    const sub = Particle.subset(list, indices);
    expect(GridPositions.is(sub.coordinates)).toBe(true);
    const grid = sub.coordinates as GridPositions;
    const source = list.coordinates as GridPositions;
    expect(grid.origin).toBe(source.origin);
    indices.forEach((src, i) => {
      expect([grid.data[i * 2], grid.data[i * 2 + 1]]).toEqual([source.data[src * 2], source.data[src * 2 + 1]]);
    });
    expect(sub.rotations).toBeInstanceOf(Int16Array);
  });
});

describe('ParticleTiles.fromList', () => {
  it('puts every particle in one tile, the one of its center, keys by tile, the tiles their targets', () => {
    const list = createList();
    const set = ParticleTiles.fromList(list, [3, 2, 2]);
    expect(set.count).toBe(Count);
    expect(set.frameCount).toBe(1);
    expect(ParticleTileGrid.cellCount(set.grid)).toBe(12);
    let total = 0;
    let base = 0;
    const seen = new Set<number>();
    const p = Vec3();
    for (const [id, tile] of set.tiles) {
      expect(tile.id).toBe(id);
      expect(tile.base).toBe(base);
      base += tile.count;
      const t = set.provider.getTile(id, 0) as ParticleList;
      expect(Task.is(t)).toBe(false);
      expect(t.count).toBe(tile.count);
      // no targets' data: the catalog's ids only
      expect(t.targetMapping).toBeUndefined();
      expect(t.targetInfo).toBe(set.catalog.targetInfo);
      expect(Array.from(t.keys!)).toEqual(Array.from({ length: tile.count }, (_, i) => tile.base + i));
      const targets = new Set<number>();
      for (let i = 0; i < t.count; ++i) {
        const src = sourceIndex(t, i);
        expect(seen.has(src)).toBe(false);
        seen.add(src);
        GridPositions.get(t.coordinates, i, p);
        expect(Box3D.containsVec3(tile.box, p)).toBe(true);
        expect(ParticleTileGrid.getIdOf(set.grid, p[0], p[1], p[2])).toBe(id);
        targets.add(ParticleValues.get(t.targets, i));
      }
      expect(Array.from(tile.targets)).toEqual([...targets].sort((a, b) => a - b));
      total += tile.count;
    }
    expect(total).toBe(Count);
    expect(seen.size).toBe(Count);
    // a tile's list once: the same object again
    const [first] = set.tiles.keys();
    expect(set.provider.getTile(first, 0)).toBe(set.provider.getTile(first, 0));
  });

  it("knows the sphere all its particles are within (its tiles' spheres'), padded by their largest radius", () => {
    const list = createList();
    const radii = ParticleValues.ofRuns(Int32Array.of(0, 50, Count), Float32Array.of(2, 4));
    for (const [l, radius] of [
      [list, 0],
      [{ ...list, radii }, 4],
    ] as const) {
      const sphere = ParticleTiles.fromList(l, [3, 2, 2]).boundingSphere!;
      const p = Vec3();
      const box = Box3D.setEmpty(Box3D());
      let farthest = 0;
      for (let i = 0; i < Count; ++i) {
        GridPositions.get(l.coordinates, i, p);
        Box3D.add(box, p);
        farthest = Math.max(farthest, Vec3.distance(p, sphere.center) + radius);
      }
      expect(farthest).toBeLessThanOrEqual(sphere.radius + 1e-6);
      // not much more: their box's sphere's at most a fifth larger
      const own = Sphere3D.fromBox3D(Sphere3D(), box);
      expect(sphere.radius).toBeLessThan((own.radius + radius) * 1.2);
    }
  });

  it("serves a tile as a view of the list in tile order: its particles' values as the list's, no copy a tile", () => {
    const list = createList();
    const set = ParticleTiles.fromList(list, [3, 2, 2]);
    const p = Vec3();
    const q = Vec3();
    const tiles = [...set.tiles.keys()].map((id) => set.provider.getTile(id, 0) as ParticleList);
    for (const t of tiles) {
      for (let i = 0; i < t.count; ++i) {
        const src = sourceIndex(t, i);
        GridPositions.get(t.coordinates, i, p);
        GridPositions.get(list.coordinates, src, q);
        expect(Vec3.distance(p, q)).toBeLessThan(1e-3);
        expect(Array.from(t.rotations!.subarray(i * 4, i * 4 + 4))).toEqual(
          Array.from(list.rotations!.subarray(src * 4, src * 4 + 4)),
        );
        expect(ParticleValues.get(t.targets, i)).toBe(ParticleValues.get(list.targets, src));
        expect(ParticleValues.get(t.entities!, i)).toBe(ParticleValues.get(list.entities!, src));
        expect(t.attributes!.get('score')!.column.value(i)).toBeCloseTo(src / 10, 5);
      }
    }
    // a tile's bounds given with it (its centers' box from the tiling's pass), holding its particles
    for (const t of tiles) {
      expect(t._propertyData[Particle.BoundaryDescriptor.name]).toBeDefined();
      const { sphere } = Particle.getBoundary(t);
      for (let i = 0; i < t.count; ++i) {
        GridPositions.get(t.coordinates, i, p);
        expect(Vec3.distance(p, sphere.center)).toBeLessThanOrEqual(sphere.radius + 1e-3);
      }
    }
    // views of one ordering: the tiles' rotations and positions share their buffers
    expect(tiles[0].rotations!.buffer).toBe(tiles[1].rotations!.buffer);
    // a grid's words, else float32 positions
    const words = (t: ParticleList) => (GridPositions.is(t.coordinates) ? t.coordinates.data : t.coordinates);
    expect(words(tiles[0]).buffer).toBe(words(tiles[1]).buffer);
  });

  it("knows the targets' kinds up front, their data served on demand as the same objects", () => {
    const list = createList();
    const set = ParticleTiles.fromList(list, [2, 2, 2]);
    expect([...set.catalog.targets.values()]).toEqual([
      {
        id: 1,
        kind: 'sites',
        sites: { quantized: false, elements: false },
        radius: 1.5,
        volume: expect.closeTo(3 * (4 / 3) * Math.PI * 1.5 ** 3 * 0.7, 6),
        byteCount: 48,
      },
      {
        id: 2,
        kind: 'sites',
        sites: { quantized: false, elements: false },
        radius: 1.5,
        volume: expect.closeTo(5 * (4 / 3) * Math.PI * 1.5 ** 3 * 0.7, 6),
        byteCount: 80,
      },
    ]);
    expect(set.catalog.entityInfo).toBe(list.entityInfo);
    expect(set.provider.getTarget(2)).toBe(list.targetMapping!.get(2));
    expect(() => set.provider.getTarget(3)).toThrow();
    // quantized with elements
    const quantized = ParticleSites.ofQuantized(1, {
      data: Int16Array.of(0, 0, 0, 0),
      origin: Vec3(),
      scale: 1,
      typeRadii: Float32Array.of(1.5),
      typeElements: Uint8Array.of(8),
    });
    const withQuantized = { ...list, targetMapping: new Map([[1, { kind: 'sites' as const, sites: quantized }]]) };
    const info = ParticleTiles.fromList(withQuantized, [1, 1, 1]).catalog.targets.get(1)!;
    expect(info.sites).toEqual({ quantized: true, elements: true });
  });

  it('serves its tiles and targets as tasks after a delay, as a fetch would', async () => {
    vi.useFakeTimers();
    try {
      const set = ParticleTiles.fromList(createList(), [2, 2, 2], { delay: 50 });
      const [id] = set.tiles.keys();
      const task = set.provider.getTile(id, 0);
      expect(Task.is(task)).toBe(true);
      const running = (task as Task<ParticleList>).run();
      await vi.advanceTimersByTimeAsync(50);
      expect((await running).count).toBe(set.tiles.get(id)!.count);
      const target = set.provider.getTarget(1);
      expect(Task.is(target)).toBe(true);
      const loading = (target as Task<ParticleTarget>).run();
      await vi.advanceTimersByTimeAsync(50);
      expect((await loading).kind).toBe('sites');
    } finally {
      vi.useRealTimers();
    }
  });
});

describe('ParticleTiles.fromTrajectory', () => {
  it("keeps a tile's particles over the frames, their poses the frame's, its box where they go", () => {
    const frames = [createList(), createList(30), createList(60)];
    const set = ParticleTiles.fromTrajectory(new ArrayParticleTrajectory(frames), [2, 2, 2]);
    expect(set.frameCount).toBe(3);
    const p = Vec3();
    const q = Vec3();
    for (const [id, tile] of set.tiles) {
      const first = set.provider.getTile(id, 0) as ParticleList;
      for (let f = 1; f < 3; ++f) {
        const t = set.provider.getTile(id, f) as ParticleList;
        // the same particles posed otherwise: one line of poses
        expect(Particle.isSamePoses(first, t)).toBe(true);
        expect(t.count).toBe(tile.count);
        for (let i = 0; i < t.count; ++i) {
          GridPositions.get(t.coordinates, i, p);
          GridPositions.get(frames[f].coordinates, sourceIndex(first, i), q);
          expect(Array.from(p)).toEqual(Array.from(q));
          expect(Box3D.containsVec3(tile.box, p)).toBe(true);
        }
      }
    }
  });
});

describe('ParticleValues.slice', () => {
  it('takes a range of runs as runs over it (cut at its ends), of an array as a view', () => {
    const runs = ParticleValues.ofRuns(Int32Array.of(0, 4, 10, 12), Int32Array.of(7, 8, 9));
    const s = ParticleValues.slice(runs, 2, 9);
    expect(ParticleValues.isRuns(s)).toBe(true);
    expect(Array.from({ length: 9 }, (_, i) => ParticleValues.get(s, i))).toEqual(
      Array.from({ length: 9 }, (_, i) => ParticleValues.get(runs, 2 + i)),
    );
    // within one run: one run
    const one = ParticleValues.slice(runs, 5, 3);
    expect(Array.from((one as { values: Int32Array }).values)).toEqual([8]);
    const array = Int32Array.of(1, 2, 3, 4, 5);
    const view = ParticleValues.slice(array, 1, 3) as Int32Array;
    expect(Array.from(view)).toEqual([2, 3, 4]);
    expect(view.buffer).toBe(array.buffer);
  });
});

describe('ParticleTileCatalog.getCompartmentEntities', () => {
  /** compartments `a` and `b` over runs of entities (-1 none), the entities' order the info's, not their ids' */
  function createCompartments(a: number, b: number): ParticleList {
    const base = createList();
    return {
      ...base,
      compartments: ParticleValues.ofRuns(Int32Array.of(0, 50, 100, 150, Count), Int32Array.of(a, b, a, b)),
      compartmentInfo: new Map([
        [a, { name: 'a' }],
        [b, { name: 'b' }],
      ]),
      entities: ParticleValues.ofRuns(
        Int32Array.of(0, 20, 50, 60, 100, 120, 150, 170, Count),
        Int32Array.of(1, 0, 1, -1, 0, 1, 1, -1),
      ),
      entityInfo: new Map([
        [1, { name: 'e1' }],
        [0, { name: 'e0' }],
      ]),
      _propertyData: Object.create(null),
    };
  }

  it("lists a compartment's entities once each, in the info's order, and whether some of its particles have none", () => {
    // ids small (a table of pairs seen) and one past it (the sets alone): the same
    for (const [a, b] of [
      [3, 7],
      [3, 1 << 23],
    ]) {
      const out = ParticleTileCatalog.getCompartmentEntities(createCompartments(a, b))!;
      expect(out.get(a)).toEqual({ name: 'a', entities: [1, 0], withoutEntity: false });
      expect(out.get(b)).toEqual({ name: 'b', entities: [1], withoutEntity: true });
    }
  });
});

describe('ParticleTiles.getAutoCells', () => {
  it("takes a tile-ordered list's written grid, else the power of two of about a page of particles a cell", () => {
    const list = createList();
    // 200 particles: one cell
    expect(ParticleTiles.getAutoCells(list)).toBe(1);
    const sized = (count: number) => ({ ...list, count, _propertyData: Object.create(null) }) as ParticleList;
    // the yeast's 31.5 M: 16 (3,845 pages' worth of cells, a cube root of 15.7)
    expect(ParticleTiles.getAutoCells(sized(31_500_000))).toBe(16);
    expect(ParticleTiles.getAutoCells(sized(8192 * 512))).toBe(8);
    expect(ParticleTiles.getAutoCells(sized(1e12))).toBe(64);
    // written tile-ordered: its grid, whatever its count
    const ordered = sized(31_500_000);
    ParticleTiles.setDirectory(ordered, {
      min: Vec3.create(0, 0, 0),
      max: Vec3.create(1, 1, 1),
      dimensions: 4,
      ranges: new Int32Array(0),
      boxes: new Float32Array(0),
    });
    expect(ParticleTiles.getAutoCells(ordered)).toBe(4);
  });
});
