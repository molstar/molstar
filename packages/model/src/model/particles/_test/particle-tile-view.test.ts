/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { CustomProperties } from '../../custom-property.js';
import {
  Particle,
  type ParticleList,
  ParticleRotations,
  type ParticleTarget,
  ParticleValues,
} from '../particle-list.js';
import { ParticleSlots } from '../particle-slots.js';
import { ParticleTileView, type ParticleTileViewSource } from '../particle-tile-view.js';
import { type ParticleTileSet, ParticleTiles } from '../particle-tiles.js';
import { ArrayParticleTrajectory } from '../particle-trajectory.js';

const Count = 120;

/** `Count` particles spread over 0-100 Å, targets 1, 2 and 3 in runs of ten, entities by target, an attribute */
function createList(shift = 0): ParticleList {
  const coordinates = new Float32Array(Count * 3);
  const rotations = new Float32Array(Count * 4);
  for (let i = 0; i < Count; ++i) {
    coordinates.set([((i * 37) % 101) + shift, (i * 53) % 97, (i * 11) % 89], i * 3);
    rotations.set([0, 0, Math.sin(i), Math.cos(i)], i * 4);
  }
  const targets = Int32Array.from({ length: Count }, (_, i) => 1 + (Math.floor(i / 10) % 3));
  const sites = (n: number): ParticleTarget => ({
    kind: 'sites',
    sites: { count: n, positions: new Float32Array(n * 3), radii: new Float32Array(n).fill(1.5) },
  });
  return {
    label: 'test',
    count: Count,
    targets,
    targetInfo: new Map([
      [1, { name: 'one' }],
      [2, { name: 'two' }],
      [3, { name: 'three' }],
    ]),
    targetMapping: new Map([
      [1, sites(2)],
      [2, sites(3)],
      [3, sites(4)],
    ]),
    entities: Int32Array.from(targets, (t) => t * 10),
    entityInfo: new Map([
      [10, { name: 'e1' }],
      [20, { name: 'e2' }],
      [30, { name: 'e3' }],
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
          max: 12,
        },
      ],
    ]),
    getParticleLabel: (index) => `particle ${index}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

/** a set's tiles and targets as a loader's ready ones */
const sourceOf = (set: ParticleTileSet): ParticleTileViewSource => ({
  getTile: (id, frame) => set.provider.getTile(id, frame) as ParticleList,
  getTarget: (id) => set.provider.getTarget(id) as ParticleTarget,
});

/** the index in the source list of a view's particle `i`, by its label (a tile's particle labelled by its index) */
const sourceIndex = (list: ParticleList, i: number) => Number(list.getParticleLabel(i).split(' ')[1]);

describe('particle tile view', () => {
  const list = createList();
  const set = ParticleTiles.fromList(list, [2, 2, 2]);
  const source = sourceOf(set);
  const ids = [...set.tiles.keys()].sort((a, b) => a - b);

  it("concatenates the tiles' particles by target, each the source's: its position, rotation, values", () => {
    const view = ParticleTileView.create(set, source, ids, 0);
    const v = view.list;
    expect(v.count).toBe(Count);
    expect(ParticleTileView.of(v)).toBe(view);
    // a target a run, in ascending order
    expect(ParticleValues.isRuns(v.targets)).toBe(true);
    const runs: number[] = [];
    ParticleValues.forEachRun(v.targets, v.count, (_s, _e, t) => runs.push(t));
    expect(runs).toEqual([1, 2, 3]);
    // every source particle once
    const seen = new Set<number>();
    const p = Vec3();
    const q = Vec3();
    for (let i = 0; i < v.count; ++i) {
      const src = sourceIndex(v, i);
      seen.add(src);
      expect(Array.from(GridPositions.get(v.coordinates, i, p))).toEqual(
        Array.from(GridPositions.get(list.coordinates, src, q)),
      );
      expect(Array.from(v.rotations!.subarray(i * 4, i * 4 + 4))).toEqual(
        Array.from(list.rotations!.subarray(src * 4, src * 4 + 4)),
      );
      expect(ParticleValues.get(v.targets, i)).toBe(ParticleValues.get(list.targets, src));
      expect(ParticleValues.get(v.entities!, i)).toBe(ParticleValues.get(list.entities!, src));
      expect(v.attributes!.get('score')!.column.value(i)).toBeCloseTo(src / 10, 5);
    }
    expect(seen.size).toBe(Count);
    // the targets the same objects, their info and the entities' the catalog's
    expect(v.targetMapping!.get(2)).toBe(list.targetMapping!.get(2));
    expect(v.targetInfo.get(2)).toBe(list.targetInfo.get(2));
    expect(v.entityInfo).toBe(list.entityInfo);
  });

  it('keys a particle by its index in the whole, the same in views of other tiles', () => {
    const a = ParticleTileView.create(set, source, ids.slice(0, 4), 0);
    const b = ParticleTileView.create(set, source, ids.slice(2), 0);
    const keyOf = (view: ParticleTileView) => {
      const out = new Map<number, number>();
      for (let i = 0; i < view.list.count; ++i) out.set(sourceIndex(view.list, i), view.list.keys![i]);
      return out;
    };
    const ka = keyOf(a);
    const kb = keyOf(b);
    let shared = 0;
    for (const [src, key] of ka) {
      if (!kb.has(src)) continue;
      shared += 1;
      expect(kb.get(src)).toBe(key);
    }
    expect(shared).toBeGreaterThan(0);
    // all keys of the whole distinct, from 0
    const all = ParticleTileView.create(set, source, ids, 0).list.keys!;
    expect([...all].sort((x, y) => x - y)).toEqual(Array.from({ length: Count }, (_, i) => i));
    // the targets the same objects in both
    expect(a.list.targetMapping!.get(1)).toBe(b.list.targetMapping!.get(1));
  });

  it("gives a target's composition: its tiles and their ranges; the same part where its tiles are", () => {
    const view = ParticleTileView.create(set, source, ids, 0);
    let end = 0;
    for (const target of [1, 2, 3]) {
      const part = view.composition.get(target)!;
      expect([...part.tiles]).toEqual([...part.tiles].sort((x, y) => x - y));
      expect(part.offsets[0]).toBe(end);
      end = part.offsets[part.tiles.length];
      for (let k = 0; k < part.tiles.length; ++k) {
        for (let i = part.offsets[k]; i < part.offsets[k + 1]; ++i) {
          expect(ParticleValues.get(view.list.targets, i)).toBe(target);
          // the particle of that tile: its key in the tile's range of the whole
          const tile = set.tiles.get(part.tiles[k])!;
          expect(view.list.keys![i]).toBeGreaterThanOrEqual(tile.base);
          expect(view.list.keys![i]).toBeLessThan(tile.base + tile.count);
        }
      }
    }
    expect(end).toBe(Count);
    const again = ParticleTileView.create(set, source, ids, 0);
    expect(ParticleTileView.isSamePart(view, again, 1)).toBe(true);
    const fewer = ParticleTileView.create(set, source, ids.slice(1), 0);
    const dropped = [1, 2, 3].filter((t) => !ParticleTileView.isSamePart(view, fewer, t));
    expect(dropped.length).toBeGreaterThan(0);
    // a target in neither: the same (none)
    expect(ParticleTileView.isSamePart(view, fewer, 99)).toBe(true);
  });

  it('steps a frame as the same particles posed otherwise', () => {
    const trajectory = new ArrayParticleTrajectory([createList(), createList(5)]);
    const tiled = ParticleTiles.fromTrajectory(trajectory, [2, 2, 2]);
    const s = sourceOf(tiled);
    const tiles = [...tiled.tiles.keys()];
    const v0 = ParticleTileView.create(tiled, s, tiles, 0);
    const v1 = ParticleTileView.create(tiled, s, tiles, 1, v0);
    expect(Particle.isSamePoses(v0.list, v1.list)).toBe(true);
    expect(v1.composition).toBe(v0.composition);
    expect(v1.frame).toBe(1);
    expect(ParticleTileView.of(v1.list)).toBe(v1);
    const p = Vec3();
    const q = Vec3();
    for (let i = 0; i < v0.list.count; ++i) {
      GridPositions.get(v0.list.coordinates, i, p);
      GridPositions.get(v1.list.coordinates, i, q);
      expect(q[0] - p[0]).toBeCloseTo(5, 4);
    }
    // without the previous view, or of other tiles: a list of its own
    expect(Particle.isSamePoses(v0.list, ParticleTileView.create(tiled, s, tiles, 1).list)).toBe(false);
    const other = ParticleTileView.create(tiled, s, tiles.slice(1), 1, v0);
    expect(Particle.isSamePoses(v0.list, other.list)).toBe(false);
  });

  it('keeps positions on one grid exact and snorm16 rotations as they are; mixed rotations as float32', () => {
    const base = createList();
    const packed = {
      ...base,
      coordinates: GridPositions.fromFloat(base.coordinates as Float32Array),
      rotations: ParticleRotations.toSnorm16(base.rotations!),
    };
    const grid = ParticleTiles.fromList(packed, [2, 2, 2]);
    const view = ParticleTileView.create(grid, sourceOf(grid), [...grid.tiles.keys()], 0);
    expect(GridPositions.is(view.list.coordinates)).toBe(true);
    expect(view.list.rotations).toBeInstanceOf(Int16Array);
    const src = packed.coordinates as GridPositions;
    const out = view.list.coordinates as GridPositions;
    for (let i = 0; i < view.list.count; ++i) {
      const k = sourceIndex(view.list, i);
      expect([out.data[i * 2], out.data[i * 2 + 1]]).toEqual([src.data[k * 2], src.data[k * 2 + 1]]);
    }
    // one tile's rotations float32, the others' snorm16
    const s = sourceOf(grid);
    const first = [...grid.tiles.keys()][0];
    const mixed: ParticleTileViewSource = {
      getTarget: s.getTarget,
      getTile: (id, f) => {
        const t = s.getTile(id, f)!;
        if (id !== first) return t;
        const r = t.rotations as Int16Array;
        return { ...t, rotations: Float32Array.from(r, (x) => x / ParticleRotations.Snorm16One) };
      },
    };
    const m = ParticleTileView.create(grid, mixed, [...grid.tiles.keys()], 0).list;
    expect(m.rotations).toBeInstanceOf(Float32Array);
    for (let i = 0; i < m.count; ++i) {
      const k = sourceIndex(m, i);
      for (let c = 0; c < 4; ++c) {
        expect(m.rotations![i * 4 + c]).toBeCloseTo(packed.rotations[k * 4 + c] / ParticleRotations.Snorm16One, 6);
      }
    }
  });

  it('labels a particle as its tile does; throws for a tile not ready; an empty view', () => {
    const view = ParticleTileView.create(set, source, ids, 0);
    expect(view.list.getParticleLabel(7)).toMatch(/^particle \d+$/);
    const notReady: ParticleTileViewSource = { getTile: () => undefined, getTarget: source.getTarget };
    expect(() => ParticleTileView.create(set, notReady, ids.slice(0, 1), 0)).toThrow(/not ready/);
    const empty = ParticleTileView.create(set, source, [], 0);
    expect(empty.list.count).toBe(0);
    expect(empty.composition.size).toBe(0);
  });
});

describe('particle tile slot view', () => {
  const list = createList();
  const set = ParticleTiles.fromList(list, [2, 2, 2]);
  const source = sourceOf(set);
  const ids = [...set.tiles.keys()].sort((a, b) => a - b);

  it("holds a tile's list as its page at its base, a particle's slot its key", () => {
    const view = ParticleTileView.createSlots(set, source, ids, 0);
    const v = view.list;
    expect(view.slots).toBe(true);
    expect(ParticleTileView.of(v)).toBe(view);
    expect(v.count).toBe(Count);
    expect(ParticleSlots.getLiveCount(v)).toBe(Count);
    for (const id of ids) {
      const tile = set.tiles.get(id)!;
      const page = v.pages!.pages[ParticleSlots.findPage(v, tile.base)];
      expect(page.start).toBe(tile.base);
      // the tile's list as it is
      expect(page.list).toBe(source.getTile(id, 0));
      for (let i = 0; i < tile.count; ++i) expect(page.list.keys![i]).toBe(tile.base + i);
    }
    // the same particles as the concatenated view's, by key
    const flat = ParticleTileView.create(set, source, ids, 0).list;
    const flattened = ParticleSlots.flatten(v);
    const p = Vec3();
    const q = Vec3();
    const byKey = new Map<number, number>();
    for (let i = 0; i < flat.count; ++i) byKey.set(flat.keys![i], i);
    for (let i = 0; i < flattened.count; ++i) {
      const j = byKey.get(flattened.keys![i])!;
      expect(Array.from(GridPositions.get(flattened.coordinates, i, p))).toEqual(
        Array.from(GridPositions.get(flat.coordinates, j, q)),
      );
    }
    expect(v.targetMapping!.get(1)).toBe(list.targetMapping!.get(1));
    // no composition: a consumer updates by its changes
    expect(ParticleTileView.isSamePart(view, view, 1)).toBe(false);
  });

  it('adds and lets go of tiles as pages, the others kept, the changes by page', () => {
    const v1 = ParticleTileView.createSlots(set, source, ids.slice(0, 3), 0);
    const v2 = ParticleTileView.createSlots(set, source, ids.slice(1, 4), 0, v1.list);
    const changes = ParticleSlots.getChanges(v1.list, v2.list)!;
    const startOf = (id: number) => set.tiles.get(id)!.base;
    expect(changes.added.map((p) => p.start)).toEqual([startOf(ids[3])]);
    expect(changes.removed.map((p) => p.start)).toEqual([startOf(ids[0])]);
    expect(changes.posed).toEqual([]);
    // the tiles in both: the same page objects
    for (const id of ids.slice(1, 3)) {
      const at = (v: typeof v1) => v.list.pages!.pages[ParticleSlots.findPage(v.list, startOf(id))];
      expect(at(v2)).toBe(at(v1));
    }
    // of another set: no line, built again
    const other = ParticleTiles.fromList(list, [2, 2, 2]);
    const v3 = ParticleTileView.createSlots(other, sourceOf(other), ids.slice(1, 4), 0, v2.list);
    expect(ParticleSlots.getChanges(v2.list, v3.list)).toBeUndefined();
  });

  it('steps a frame as its pages posed otherwise', () => {
    const trajectory = new ArrayParticleTrajectory([createList(), createList(5)]);
    const tiled = ParticleTiles.fromTrajectory(trajectory, [2, 2, 2]);
    const s = sourceOf(tiled);
    const tiles = [...tiled.tiles.keys()];
    const v0 = ParticleTileView.createSlots(tiled, s, tiles, 0);
    const v1 = ParticleTileView.createSlots(tiled, s, tiles, 1, v0.list);
    const changes = ParticleSlots.getChanges(v0.list, v1.list)!;
    expect(changes.added).toEqual([]);
    expect(changes.removed).toEqual([]);
    expect(changes.posed.length).toBe(tiles.length);
  });
});
