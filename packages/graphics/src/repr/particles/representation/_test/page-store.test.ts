/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { Box3D } from '@molstar/core/math/geometry/primitives/box3d';
import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { CustomProperties } from '@molstar/model/model/custom-property';
import {
  Particle,
  type ParticleList,
  ParticleRotations,
  type ParticleTarget,
  ParticleValues,
} from '@molstar/model/model/particles/particle-list';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import { ParticleTiles } from '@molstar/model/model/particles/particle-tiles';
import {
  CopyPoseFloats,
  getStoreFormatOf,
  ExtraWords,
  PageFlagsWord,
  TemplatePage,
  writeCopyPose,
  NoMember,
  PageRecordWords,
  ParticlePageStore,
  PerParticleMembers,
  StorePageSize,
  type StorePlacement,
} from '../target/page-store.js';

/** `n` particles, targets by `targetOf`, rotated by `i` about an axis of their own, radii a particle or by run */
function createList(n: number, targetOf: (i: number) => number, radii?: 'particle' | 'runs'): ParticleList {
  const coordinates = new Float32Array(n * 3);
  const rotations = new Float32Array(n * 4);
  const q = Quat();
  for (let i = 0; i < n; ++i) {
    coordinates.set([i * 0.5, 2 * i - 7, 3 + (i % 5)], i * 3);
    Quat.setAxisAngle(q, Vec3.normalize(Vec3(), Vec3.create(1, i % 3, 2)), i * 0.3);
    rotations.set([q[0], q[1], q[2], q[3]], i * 4);
  }
  const targets = Int32Array.from({ length: n }, (_, i) => targetOf(i));
  return {
    count: n,
    targets,
    targetInfo: new Map(),
    coordinates,
    rotations,
    radii:
      radii === 'particle'
        ? Float32Array.from({ length: n }, (_, i) => 1 + (i % 7) * 0.25)
        : radii === 'runs'
          ? ParticleValues.ofRuns(Int32Array.of(0), Float32Array.of(1.5))
          : undefined,
    getParticleLabel: (i) => `${i}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

const centers = new Map<number, Vec3>([
  [1, Vec3.create(1, 0, 0)],
  [2, Vec3.create(0, 2, 0)],
  [3, Vec3.create(0, 0, -1)],
  [4, Vec3.create(0.5, 0.5, 0.5)],
]);
const placement: StorePlacement = {
  centerOf: (t) => centers.get(t)!,
  // target 2 scaled by its radii; target 4 not the store's
  scaleByRadius: (t) => t === 2,
  accepts: (t) => t !== 4,
};

/** floats of a composed pose: `t` and the scale, the quaternion */
const PoseFloats = 8;

/**
 * the instances a store makes of `list` from slot `start`, by slot: a particle of a target `placement` accepts, its
 * pose `t = p - s R c` (`c` its target's center, `s` its radius where the target scales by radius, else 1), `s`,
 * its quaternion
 */
function expected(list: ParticleList, start: number) {
  const out = new Map<number, { target: number; pose: Float32Array }>();
  const coordinates = list.coordinates as Float32Array;
  const p = Vec3();
  const rc = Vec3();
  for (let i = 0; i < list.count; ++i) {
    const target = list.targets[i];
    if (placement.accepts && !placement.accepts(target)) continue;
    Vec3.set(p, coordinates[i * 3], coordinates[i * 3 + 1], coordinates[i * 3 + 2]);
    const r = list.rotations as Float32Array;
    const q = Quat.create(r[i * 4], r[i * 4 + 1], r[i * 4 + 2], r[i * 4 + 3]);
    const scale =
      (placement.scaleByRadius?.(target) && list.radii ? ParticleValues.get(list.radii, i) : 1) *
      (list.scales ? ParticleValues.get(list.scales, i) : (list.scale ?? 1));
    Vec3.transformQuat(rc, placement.centerOf(target), q);
    Vec3.scaleAndAdd(p, p, rc, -scale);
    out.set(start + i, { target, pose: Float32Array.of(p[0], p[1], p[2], scale, q[0], q[1], q[2], q[3]) });
  }
  return out;
}

/** `store`'s instances as `want`'s: by slot, its member its target's, its pose within `tolerance` */
function expectInstances(
  store: ParticlePageStore,
  keys: number[],
  want: Map<number, { target: number; pose: Float32Array }>,
  tolerance = 1e-4,
) {
  let n = 0;
  for (const key of keys) {
    const { slots, poses, members } = store.compose(store.getBlock(key)!);
    for (let i = 0; i < slots.length; ++i) {
      const w = want.get(slots[i])!;
      expect(w, `slot ${slots[i]}`).toBeDefined();
      expect(members[i]).toBe(store.pool.get(w.target));
      for (let k = 0; k < PoseFloats; ++k)
        expect(Math.abs(poses[i * PoseFloats + k] - w.pose[k])).toBeLessThan(tolerance);
      n += 1;
    }
  }
  expect(n).toBe(want.size);
}

describe('particle page store', () => {
  it("holds a list's particles in pages: a target a member, its count its instances', the poses by slot", () => {
    for (const radii of ['particle', 'runs', undefined] as const) {
      const pages = new ParticlePageStore();
      const a = createList(40, (i) => (i < 10 ? 1 : i < 18 ? 4 : i < 30 ? 2 : 3), radii);
      const b = createList(25, (i) => (i % 2 ? 2 : 3), radii);
      pages.add(1, a, 100, placement);
      pages.add(2, b, 200, placement);
      const want = new Map([...expected(a, 100), ...expected(b, 200)]);
      // target 4 not the store's: no member, its particles none
      expect(pages.pool.get(4)).toBe(-1);
      for (const t of [1, 2, 3]) {
        const count = [...want.values()].filter((w) => w.target === t).length;
        expect(pages.pool.counts[pages.pool.get(t)]).toBe(count);
      }
      expectInstances(pages, [1, 2], want);
    }
  });

  it("draws a list's particles larger by its scale (a coarse level's node), the store's largest scale kept", () => {
    for (const radii of ['particle', 'runs', undefined] as const) {
      const pages = new ParticlePageStore();
      const a = { ...createList(40, (i) => (i < 10 ? 1 : i < 30 ? 2 : 3), radii), scale: 4 };
      const b = createList(25, (i) => (i % 2 ? 2 : 3), radii);
      pages.add(1, a, 100, placement);
      expect(pages.maxScale).toBe(4);
      pages.add(2, b, 200, placement);
      expect(pages.maxScale).toBe(4);
      expectInstances(pages, [1, 2], new Map([...expected(a, 100), ...expected(b, 200)]));
    }
  });

  it("draws a list's particles larger by their runs' scales (a node by target), in place of its scale", () => {
    for (const radii of ['particle', 'runs', undefined] as const) {
      const pages = new ParticlePageStore();
      const scales = ParticleValues.ofRuns(Int32Array.of(0, 10, 30, 40), Float32Array.of(1.5, 3, 2));
      const a = { ...createList(40, (i) => (i < 10 ? 1 : i < 30 ? 2 : 3), radii), scale: 4, scales };
      pages.add(1, a, 100, placement);
      // the largest run's, not the list's
      expect(pages.maxScale).toBe(3);
      expectInstances(pages, [1], expected(a, 100));
    }
  });

  it('splits a list larger than a page: its pages its slots from its start, its runs at their ends', () => {
    const store = new ParticlePageStore();
    const n = StorePageSize + 1000;
    // a run across the pages' end
    const list = createList(n, (i) => (i < StorePageSize - 10 ? 1 : 3));
    const block = store.add(7, list, 5000, placement);
    expect(block.pages.length).toBe(2);
    const [p0, p1] = block.pages;
    const r = store.records;
    expect(r[p0 * PageRecordWords + 6]).toBe(StorePageSize);
    expect(r[p0 * PageRecordWords + 7]).toBe(5000);
    expect(r[p1 * PageRecordWords + 6]).toBe(1000);
    expect(r[p1 * PageRecordWords + 7]).toBe(5000 + StorePageSize);
    // the first page two runs (the second cut at its end), the second one
    expect(r[p0 * PageRecordWords + 5]).toBe(2);
    expect(r[p1 * PageRecordWords + 5]).toBe(1);
    const first = r[p1 * PageRecordWords + 4];
    expect(Array.from(store.runs.subarray(first * 4, first * 4 + 3))).toEqual([0, 1000, store.pool.get(3)]);
    // its box: the positions' extent
    const f = store.recordFloats;
    expect(f[p1 * PageRecordWords + 8]).toBeCloseTo(StorePageSize * 0.5, 2);
    expect(f[p1 * PageRecordWords + 12]).toBeCloseTo((n - 1) * 0.5, 2);
    // every particle an instance, by its slot
    const { slots } = store.compose(block);
    expect(slots.length).toBe(n);
    expect(slots[n - 1]).toBe(5000 + n - 1);
  });

  it('keeps a page of many runs as a member a particle (its runs more than an eighth of its particles)', () => {
    const pages = new ParticlePageStore();
    // targets alternating: a run a particle
    const list = createList(64, (i) => [1, 2, 4][i % 3], 'particle');
    const block = pages.add(1, list, 0, placement);
    const p = block.pages[0];
    expect(pages.records[p * PageRecordWords + 5]).toBe(PerParticleMembers);
    // its extras: a member and a scale a particle, target 4's particles no member
    const e = pages.records[p * PageRecordWords + 11];
    expect(pages.extras[(e + 2) * ExtraWords]).toBe(NoMember);
    expect(pages.extras[e * ExtraWords]).toBe(pages.pool.get(1));
    expectInstances(pages, [1], expected(list, 0));
  });

  it('lets go of a block: its pages and runs free (reused), its members freed with their last instance', () => {
    const store = new ParticlePageStore({ pageCapacity: 4 });
    store.add(
      1,
      createList(10, () => 1),
      0,
      placement,
    );
    const second = store.add(
      2,
      createList(40, (i) => (i < 20 ? 1 : 2)),
      10,
      placement,
    );
    const changes = store.takeChanges();
    expect(changes.written.length).toBe(2);
    expect(changes.runs.length).toBe(4);
    const m2 = store.pool.get(2);
    store.remove(2);
    const removed = store.takeChanges();
    expect(removed.released).toEqual(Array.from(second.pages));
    expect(removed.freed).toEqual([m2]);
    expect(store.pool.get(2)).toBe(-1);
    expect(store.pool.counts[store.pool.get(1)]).toBe(10);
    // its record cleared, its page reused
    expect(store.records[second.pages[0] * PageRecordWords + 6]).toBe(0);
    const third = store.add(
      3,
      createList(4, () => 2),
      20,
      placement,
    );
    expect(third.pages[0]).toBe(second.pages[0]);
    // target 2 its member again (the one it had)
    expect(store.pool.get(2)).toBe(m2);
  });

  it("stores its first list's format: a grid's words and snorm16 rotations, a later float list converted", () => {
    const pages = new ParticlePageStore();
    const float = createList(30, (i) => (i < 15 ? 1 : 3));
    const grid: ParticleList = {
      ...float,
      coordinates: GridPositions.fromFloat(float.coordinates as Float32Array),
      rotations: Int16Array.from(float.rotations!, (x) => Math.round(x * 32767)),
    };
    pages.add(1, grid, 0, placement);
    pages.add(2, float, 100, placement);
    expect(pages.format).toEqual({ positions: 'grid', rotations: 'snorm16' });
    const p = pages.getBlock(2)!.pages[0];
    expect(pages.values[p]!.positions).toBeInstanceOf(Uint32Array);
    expect(pages.values[p]!.rotations).toBeInstanceOf(Int16Array);
    // a grid's step and snorm16's precision
    expectInstances(pages, [1, 2], new Map([...expected(float, 0), ...expected(float, 100)]), 2e-3);
  });

  it('grows its pages and runs, what is there kept', () => {
    const store = new ParticlePageStore({ pageCapacity: 4 });
    for (let k = 0; k < 6; ++k) {
      store.add(
        k,
        createList(300, (i) => 1 + (Math.floor(i / 3) % 3)),
        k * 300,
        placement,
      );
    }
    expect(store.takeChanges().grown).toBe(true);
    expect(store.pageCapacity).toBeGreaterThanOrEqual(6);
    for (let k = 0; k < 6; ++k) {
      const { slots } = store.compose(store.getBlock(k)!);
      expect(slots.length).toBe(300);
      expect(slots[0]).toBe(k * 300);
    }
  });

  it("gives a target of a large id (past the members' table) its member as a small one's, again after it is freed", () => {
    const large = (1 << 23) + 1;
    const far = (1 << 22) - 1;
    const of = new Map<number, Vec3>([
      [1, Vec3.create(1, 0, 0)],
      [far, Vec3.create(0, 1, 0)],
      [large, Vec3.create(0, 0, 1)],
    ]);
    const own: StorePlacement = { centerOf: (t) => of.get(t)! };
    const store = new ParticlePageStore();
    const list = createList(30, (i) => (i < 10 ? 1 : i < 20 ? large : far));
    store.add(1, list, 0, own);
    const members = [1, large, far].map((t) => store.pool.get(t));
    expect(new Set(members).size).toBe(3);
    for (const m of members) expect(store.pool.counts[m]).toBe(10);
    expect([...store.getBlock(1)!.members]).toEqual(members.flatMap((m) => [m, 10]));
    store.remove(1);
    for (const t of [1, large, far]) expect(store.pool.get(t)).toBe(-1);
    // back: each its member as before (its values its own)
    store.add(2, list, 0, own);
    expect([1, large, far].map((t) => store.pool.get(t))).toEqual(members);
  });

  it("takes a box given for every page of a block (a tile's), else its pages' own of their particles", () => {
    const store = new ParticlePageStore();
    const n = StorePageSize + 500;
    const given = Box3D.create(Vec3.create(-10, -20, -30), Vec3.create(4000, 3000, 2000));
    const a = store.add(
      1,
      createList(n, () => 1),
      0,
      placement,
      given,
    );
    const f = store.recordFloats;
    expect(a.pages.length).toBe(2);
    for (const p of a.pages) {
      const o = p * PageRecordWords;
      expect([f[o + 8], f[o + 9], f[o + 10], f[o + 12], f[o + 13], f[o + 14]]).toEqual([
        -10, -20, -30, 4000, 3000, 2000,
      ]);
    }
    // none given: the particles' (the first page's x from 0, its last x half its last index)
    const b = store.add(
      2,
      createList(n, () => 1),
      n,
      placement,
    );
    const o = b.pages[0] * PageRecordWords;
    expect(f[o + 8]).toBeCloseTo(0, 2);
    expect(f[o + 12]).toBeCloseTo((StorePageSize - 1) * 0.5, 2);
  });
});

describe('particle page store of a nested list', () => {
  /** three copies of a template of ten particles (targets 1 and 3), one turned, and a leaf (target 2) */
  function createNested() {
    const template = createList(10, (i) => (i < 6 ? 1 : 3));
    const top = createList(4, (i) => (i === 2 ? 2 : 0));
    const turn = Quat.setAxisAngle(Quat(), Vec3.create(0, 0, 1), Math.PI / 2);
    (top.rotations as Float32Array).set([0, 0, 0, 1, ...turn, 0, 0, 0, 1, 0, 0, 0, 1]);
    (top.targetInfo as Map<number, object>).set(0, { name: 'copy' }).set(2, { name: 'leaf' });
    top.targetMapping = new Map<number, ParticleTarget>([
      [0, { kind: 'particles', particles: template, center: Vec3.create(1, 2, 3) }],
      [2, { kind: 'shape', shape: { name: 'leaf' } as never }],
    ]);
    return { template, top };
  }

  it("holds a template's particles once, its pages its copies' alone, its copies' poses its instances' composed", () => {
    const { template, top } = createNested();
    const ids = ParticleNesting.getIdSpace(top);
    const base = ids.moves.get(template)!.target;
    // the leaf's target 2, the template's past the top list's ids: 3 + 1 and 3 + 3
    expect(base).toBe(3);
    const flatCenters = new Map<number, Vec3>([
      [2, Vec3.create(0, 0, 0)],
      [base + 1, Vec3.create(1, 0, 0)],
      [base + 3, Vec3.create(0, 0, -1)],
    ]);
    const flatPlacement: StorePlacement = { centerOf: (t) => flatCenters.get(t)!, accepts: (t) => flatCenters.has(t) };
    const store = new ParticlePageStore();
    store.add(0, top, 0, flatPlacement);
    const leaves = ParticleNesting.getLeafList(top, template);
    const block = store.add(-1, leaves, 0, flatPlacement, undefined, true);
    expect(block.template).toBe(true);
    // the top list's page not a template's, the template's one
    expect(store.records[store.getBlock(0)!.pages[0] * PageRecordWords + PageFlagsWord]).toBe(0);
    expect(store.records[block.pages[0] * PageRecordWords + PageFlagsWord]).toBe(TemplatePage);
    const tops = Int32Array.from(ParticleNesting.getCopiesOf(top, template));
    expect(Array.from(tops)).toEqual([0, 1, 3]);
    const poses = new Float32Array(tops.length * CopyPoseFloats);
    const p = Vec3();
    const q = Quat();
    for (let c = 0; c < tops.length; ++c) {
      GridPositions.get(top.coordinates, tops[c], p);
      ParticleRotations.get(top.rotations, tops[c], q);
      writeCopyPose(poses, c * CopyPoseFloats, p, q, Vec3.create(1, 2, 3));
    }
    // the template's members counted once as added, then once a copy
    const m1 = store.pool.get(base + 1);
    const m3 = store.pool.get(base + 3);
    expect([store.pool.counts[m1], store.pool.counts[m3]]).toEqual([6, 4]);
    store.setCopies(-1, { tops, poses });
    expect([store.pool.counts[m1], store.pool.counts[m3]]).toEqual([18, 12]);
    expect(store.getCopyFactor(-1)).toBe(3);
    expect(store.takeChanges().copies).toEqual([-1]);
    // every copy's instance of a template particle at the flattened list's position and rotation
    const flat = ParticleNesting.flatten(top);
    const composed = store.composeCopies(-1);
    expect(composed.tops.length).toBe(3 * 10);
    const r = Quat();
    const at = Vec3();
    for (let k = 0; k < composed.tops.length; ++k) {
      const i = ParticleNesting.getExpandedIndex(top, composed.tops[k], composed.slots[k]);
      expect(ParticleValues.get(flat.targets, i)).toBe(store.pool.targets[composed.members[k]]);
      // the pose places its member's center at the instance: t + s R c
      const o = k * 8;
      const m = composed.members[k];
      const c = Vec3.create(store.pool.centers[m * 3], store.pool.centers[m * 3 + 1], store.pool.centers[m * 3 + 2]);
      Quat.set(r, composed.poses[o + 4], composed.poses[o + 5], composed.poses[o + 6], composed.poses[o + 7]);
      Vec3.add(
        at,
        Vec3.create(composed.poses[o], composed.poses[o + 1], composed.poses[o + 2]),
        Vec3.scale(c, Vec3.transformQuat(c, c, r), composed.poses[o + 3]),
      );
      expect(Vec3.distance(at, GridPositions.get(flat.coordinates, i, Vec3()))).toBeLessThan(1e-4);
      expect(Math.abs(Quat.dot(r, ParticleRotations.get(flat.rotations, i, Quat())))).toBeCloseTo(1, 5);
    }
    // the copies let go with their template
    store.remove(-1);
    expect(store.getCopies(-1)).toBeUndefined();
    const changes = store.takeChanges();
    expect(changes.copies).toEqual([-1]);
    // its members' every instance gone: freed
    expect(changes.freed.sort()).toEqual([m1, m3].sort());
    expect(store.pool.get(base + 1)).toBe(-1);
  });

  it("takes the format of all its lists: a top list without rotations before a template's of them", () => {
    const top = { ...createList(4, () => 2), rotations: undefined };
    const template = createList(40, () => 1);
    const grid = { ...template, coordinates: GridPositions.fromFloat(template.coordinates as Float32Array) };
    const format = getStoreFormatOf([top, grid]);
    // the largest's positions, the rotated one's rotations
    expect(format).toEqual({ positions: 'grid', rotations: 'float32' });
    const store = new ParticlePageStore();
    store.format = format;
    store.add(0, top, 0, { centerOf: () => Vec3() });
    store.add(-1, grid, 0, { centerOf: () => Vec3() }, undefined, true);
    // the top list's identities
    expect(Array.from(store.values[store.getBlock(0)!.pages[0]]!.rotations!.subarray(0, 4))).toEqual([0, 0, 0, 1]);
    // its own format (the first list's) refuses the template's rotations
    const first = new ParticlePageStore();
    first.add(0, top, 0, { centerOf: () => Vec3() });
    expect(() => first.add(-1, grid, 0, { centerOf: () => Vec3() }, undefined, true)).toThrow(/rotations/);
  });

  it('refuses copies of a block not a template, a pose missing', () => {
    const { top } = createNested();
    const store = new ParticlePageStore();
    store.add(0, top, 0, { centerOf: () => Vec3(), accepts: (t) => t === 2 });
    expect(() => store.setCopies(0, { tops: Int32Array.of(0), poses: new Float32Array(8) })).toThrow(/template/);
    store.add(
      -1,
      createList(4, () => 1),
      0,
      { centerOf: () => Vec3() },
      undefined,
      true,
    );
    expect(() => store.setCopies(-1, { tops: Int32Array.of(0, 1), poses: new Float32Array(8) })).toThrow(/a copy pose/);
    expect(() => store.setCopies(-1, { tops: new Int32Array(0), poses: new Float32Array(0) })).toThrow(/no copies/);
  });

  it("reads a deferred list's page when asked, keeps its values given up within the cache's bytes", async () => {
    // three pages of particles on a grid, no rotations, deferred: read by range (counted), a tile a page
    const n = 3 * StorePageSize;
    const xyz = new Float32Array(n * 3);
    for (let i = 0; i < n; ++i) xyz.set([Math.floor(i / StorePageSize) * 100 + (i % 50), i % 7, i % 11], i * 3);
    const full = GridPositions.fromFloat(xyz);
    const list = {
      ...createList(n, () => 1),
      coordinates: { ...full, data: new Uint32Array(0) },
      rotations: undefined,
    };
    const reads: number[] = [];
    Particle.setDeferred(list, {
      read: async (first, count) => {
        reads.push(first);
        return { positions: full.data.slice(first * 2, (first + count) * 2) };
      },
    });
    ParticleTiles.setDirectory(list, {
      min: Vec3.create(0, 0, 0),
      max: Vec3.create(300, 10, 10),
      dimensions: 4,
      // a tile's cell, first particle, count; its box
      ranges: Int32Array.of(0, 0, StorePageSize, 1, StorePageSize, StorePageSize, 2, 2 * StorePageSize, StorePageSize),
      boxes: Float32Array.of(0, 0, 0, 49, 6, 10, 100, 0, 0, 149, 6, 10, 200, 0, 0, 249, 6, 10),
    });
    const pageBytes = StorePageSize * 8;
    const store = new ParticlePageStore({ valueCacheBytes: pageBytes });
    const { pages } = store.add(1, list, 0, placement);
    expect(pages.length).toBe(3);
    expect([...pages].map((p) => store.values[p])).toEqual([undefined, undefined, undefined]);
    // a page's box its tile's
    const f = store.recordFloats;
    const o = pages[1] * PageRecordWords;
    expect([f[o + 8], f[o + 12]]).toEqual([100, 149]);

    // read once while in flight, its values the grid's words
    const [a, b] = [pages[0], pages[1]];
    const reading = store.requestValues(a);
    expect(store.requestValues(a)).toBe(reading);
    await reading;
    expect(store.values[a]!.positions).toEqual(full.data.subarray(0, StorePageSize * 2));
    expect(store.requestValues(a)).toBeUndefined();
    await store.requestValues(b);
    expect(reads).toEqual([0, StorePageSize]);
    expect(store.readValueBytes).toBe(2 * pageBytes);

    // given up: kept, set again at once (no read)
    store.releaseValues(a);
    expect(store.values[a]).toBeUndefined();
    expect(store.cachedValueBytes).toBe(pageBytes);
    await store.requestValues(a);
    expect(store.values[a]).toBeDefined();
    expect(store.cachedValueBytes).toBe(0);
    expect(reads.length).toBe(2);

    // past the cache's bytes the least recently given up let go: read again
    store.releaseValues(a);
    store.releaseValues(b);
    expect(store.cachedValueBytes).toBe(pageBytes);
    await store.requestValues(b);
    expect(reads.length).toBe(2);
    await store.requestValues(a);
    expect(reads).toEqual([0, StorePageSize, 0]);

    // a list's values not deferred not given up; a block removed its kept values let go
    store.releaseValues(a);
    store.remove(1);
    expect(store.cachedValueBytes).toBe(0);
    const other = new ParticlePageStore();
    const [q] = other.add(
      2,
      createList(10, () => 1),
      0,
      placement,
    ).pages;
    other.releaseValues(q);
    expect(other.values[q]).toBeDefined();
  });
});
