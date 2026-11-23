/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { Task } from '@molstar/core/task';
import { CustomProperties } from '../../../model/custom-property.js';
import {
  Particle,
  type ParticleList,
  type ParticleSites,
  type ParticleTarget,
  ParticleValues,
} from '../../../model/particles/particle-list.js';
import { ParticleTileCatalog, ParticleTiles } from '../../../model/particles/particle-tiles.js';
import { createParticleListFromPack, readParticlePack, writeParticlePack } from '../pack.js';
import {
  createParticleTileSetFromPack,
  isParticlePackListStreamable,
  isParticlePackStreamable,
  openParticleListFromPackStream,
  openParticlePackStream,
  ParticlePackSource,
} from '../pack-stream.js';
import { ParticleNesting } from '../../../model/particles/particle-nesting.js';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';

/** 300 particles over a box, three targets of sites in an alternating order, entities, compartments, an attribute */
function createList(): ParticleList {
  const n = 300;
  const coordinates = new Float32Array(n * 3);
  for (let i = 0; i < n; ++i) coordinates.set([(i * 37) % 101, (i * 53) % 97, (i * 11) % 89], i * 3);
  const targets = Int32Array.from({ length: n }, (_, i) => [3, 1, 2][i % 3]);
  const sites = (count: number, r: number): ParticleTarget => ({
    kind: 'sites',
    sites: {
      count,
      positions: Float32Array.from({ length: count * 3 }, (_, k) => (k % 3) + count),
      radii: new Float32Array(count).fill(r),
    },
  });
  return {
    label: 'spread',
    count: n,
    targets,
    targetInfo: new Map([
      [1, { name: 'one', entity: 0 }],
      [2, { name: 'two', entity: 1 }],
      [3, { name: 'three', entity: -1 }],
    ]),
    targetMapping: new Map([
      [1, sites(4, 3)],
      [2, sites(2, 5)],
      [3, sites(1, 8)],
    ]),
    entities: Int32Array.from(targets, (t) => (t === 1 ? 0 : t === 2 ? 1 : -1)),
    entityInfo: new Map([
      [0, { name: 'e0' }],
      [1, { name: 'e1' }],
    ]),
    compartments: Int32Array.from(targets, (t) => (t === 2 ? 1 : 0)),
    compartmentInfo: new Map([
      [0, { name: 'c0' }],
      [1, { name: 'c1' }],
    ]),
    coordinates,
    rotations: new Float32Array(n * 4).map((_, k) => (k % 4 === 3 ? 1 : 0)),
    radii: Float32Array.from(targets, (t) => t * 2),
    attributes: new Map([
      [
        'source',
        {
          column: Column.ofFloatArray(Float32Array.from({ length: n }, (_, i) => i)),
          label: 'Source',
          min: 0,
          max: n - 1,
        },
      ],
    ]),
    getParticleLabel: (i) => `p${i}`,
    sourceData: { kind: 'test', name: 'test', data: undefined },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

const resolve = async <T>(v: T | Task<T>) => (Task.is(v) ? await v.run() : v);

function sitesOf(target: ParticleTarget) {
  if (target.kind !== 'sites') throw new Error('not sites');
  return target.sites as ParticleSites;
}

describe('particle pack read by range', () => {
  it("makes a tile-ordered pack's tile set as the whole read makes it, its tiles and targets read when asked", async () => {
    const bytes = writeParticlePack(createList(), { tiles: 4 });
    const whole = createParticleListFromPack(readParticlePack(bytes));
    const expected = ParticleTiles.fromList(whole, [4, 4, 4]);

    const source = ParticlePackSource.fromBytes(bytes);
    const read = vi.spyOn(source, 'read');
    const stream = await openParticlePackStream(source);
    expect(isParticlePackStreamable(stream.document)).toBe(true);
    const set = await createParticleTileSetFromPack(stream);
    // the directory, the opacities and the runs read at once, no tile's particles nor a target's sites
    const readsAtOnce = read.mock.calls.length;
    const blocksOf = (names: string[]) =>
      stream.document.blocks.filter((b) => names.includes(b.name)).map((b) => b.byteOffset);
    for (const offset of blocksOf(['position', 'rotation', 'site'])) {
      expect(read.mock.calls.some(([o]) => o === offset)).toBe(false);
    }

    expect(set.count).toBe(expected.count);
    expect(set.grid).toEqual(expected.grid);
    expect([...set.tiles.keys()]).toEqual([...expected.tiles.keys()]);
    for (const [id, tile] of expected.tiles) {
      const t = set.tiles.get(id)!;
      expect([t.count, t.base, Array.from(t.targets)]).toEqual([tile.count, tile.base, Array.from(tile.targets)]);
      expect(t.box).toEqual(tile.box);
      expect(Array.from(t.opacity!.density)).toEqual(Array.from(tile.opacity!.density));
    }
    expect(set.boundingSphere!.radius).toBeCloseTo(expected.boundingSphere!.radius, 3);
    // the catalog: the targets' reach and volume as written, the compartments' entities as from the particles
    for (const [id, info] of expected.catalog.targets) {
      const t = set.catalog.targets.get(id)!;
      expect(t.radius).toBeCloseTo(info.radius!, 2);
      expect(t.volume).toBeCloseTo(info.volume!, 2);
    }
    expect(set.catalog.compartmentInfo).toEqual(ParticleTileCatalog.getCompartmentEntities(whole));

    // a tile's list: its range's particles as the whole's tile has them (its positions' words, its runs' values)
    for (const id of expected.tiles.keys()) {
      const a = await resolve(expected.provider.getTile(id, 0));
      const b = await resolve(set.provider.getTile(id, 0));
      expect(b.count).toBe(a.count);
      expect(Array.from(b.keys!)).toEqual(Array.from(a.keys!));
      expect(GridPositions.is(b.coordinates)).toBe(true);
      expect(Array.from((b.coordinates as GridPositions).data)).toEqual(
        Array.from((a.coordinates as GridPositions).data),
      );
      expect(Array.from(b.rotations!)).toEqual(Array.from(a.rotations!));
      for (const k of ['targets', 'entities', 'compartments', 'radii'] as const) {
        expect(Array.from(ParticleValues.toArray(b[k]!, b.count, Float64Array))).toEqual(
          Array.from(ParticleValues.toArray(a[k]!, a.count, Float64Array)),
        );
      }
      const column = b.attributes!.get('source')!.column;
      expect(Array.from({ length: b.count }, (_, i) => column.value(i))).toEqual(
        Array.from({ length: a.count }, (_, i) => a.attributes!.get('source')!.column.value(i)),
      );
      expect(Particle.getBoundary(b).sphere.radius).toBeCloseTo(Particle.getBoundary(a).sphere.radius, 3);
      expect(b.getParticleLabel(0)).toBe(a.getParticleLabel(0));
    }
    expect(read.mock.calls.length).toBeGreaterThan(readsAtOnce);

    // a target's sites read once, its quantized words the whole's
    for (const id of [1, 2, 3]) {
      const a = sitesOf(await resolve(expected.provider.getTarget(id)));
      const before = read.mock.calls.length;
      const b = sitesOf(await resolve(set.provider.getTarget(id)));
      await resolve(set.provider.getTarget(id));
      expect(read.mock.calls.length).toBe(before + 1);
      expect(Array.from(b.quantized!.data)).toEqual(Array.from(a.quantized!.data));
      expect(b.count).toBe(a.count);
    }
  });

  it("keeps every siteStride-th of a target's sites (siteStride), its catalog's bytes the kept ones'", async () => {
    const bytes = writeParticlePack(createList(), { tiles: 4 });
    const whole = await createParticleTileSetFromPack(
      await openParticlePackStream(ParticlePackSource.fromBytes(bytes)),
    );
    const set = await createParticleTileSetFromPack(await openParticlePackStream(ParticlePackSource.fromBytes(bytes)), {
      siteStride: 2,
    });
    for (const id of [1, 2, 3]) {
      const a = sitesOf(await resolve(whole.provider.getTarget(id)));
      const b = sitesOf(await resolve(set.provider.getTarget(id)));
      expect(b.stride).toBe(2);
      expect(b.count).toBe(Math.ceil(a.count / 2));
      const words = Array.from({ length: b.count }, (_, i) => Array.from(a.quantized!.data.subarray(i * 8, i * 8 + 4)));
      expect(
        Array.from({ length: b.count }, (_, i) => Array.from(b.quantized!.data.subarray(i * 4, i * 4 + 4))),
      ).toEqual(words);
      expect(set.catalog.targets.get(id)!.byteCount).toBe(b.count * 8);
    }
  });

  it("writes coarse levels, an octree of samples of its tiles' particles, read as tiles of the set", async () => {
    const list = createList();
    const bytes = writeParticlePack(list, { tiles: 4, levels: true });
    // the whole read: the particles alone (the nodes' copies not read)
    const whole = createParticleListFromPack(readParticlePack(bytes));
    expect(whole.count).toBe(list.count);
    const stream = await openParticlePackStream(ParticlePackSource.fromBytes(bytes));
    expect(stream.document.tiles!.levels).toBe(2);
    const set = await createParticleTileSetFromPack(stream);
    const levels = set.levels!;
    expect(levels.count).toBe(2);
    const base = levels.base;
    expect(base.count).toBe(list.count);
    // a particle's tile, by its index in the pack
    const tileOf = new Int32Array(list.count);
    for (const t of base.tiles.values()) tileOf.fill(t.id, t.base, t.base + t.count);
    const words = (l: ParticleList, i: number) => {
      const data = (l.coordinates as GridPositions).data;
      return [data[i * 2], data[i * 2 + 1]];
    };
    const slots = new Set<number>();
    for (let level = 1; level <= 2; ++level) {
      const d = 4 >> level;
      for (let cell = 0; cell < d * d * d; ++cell) {
        const id = levels.getNode(level, cell);
        // its tiles: those in its cell of the level's grid; their particles by target
        let n = 0;
        const byTarget = new Map<number, number>();
        for (const t of base.tiles.values()) {
          if (ParticleTiles.getLevelCell(4, t.id, level) !== cell) continue;
          n += t.count;
          for (let i = t.base; i < t.base + t.count; ++i) {
            const target = ParticleValues.get(whole.targets, i);
            byTarget.set(target, (byTarget.get(target) ?? 0) + 1);
          }
        }
        if (n === 0) {
          expect(id).toBeUndefined();
          continue;
        }
        const node = set.tiles.get(id!)!;
        // of each target every 8^level-th of its particles there (at least one), its slots after the particles'
        let kept = 0;
        for (const c of byTarget.values()) kept += Math.ceil(c / 8 ** level);
        expect(node.count).toBe(kept);
        expect(node.base).toBeGreaterThanOrEqual(base.count);
        const nodeList = await resolve(set.provider.getTile(id!, 0));
        expect(nodeList.count).toBe(node.count);
        // each target's sample drawn larger by its own scale: its volume kept (its particles over those kept, cubed)
        expect(nodeList.scale).toBeUndefined();
        const keptOf = new Map<number, number>();
        for (let i = 0; i < nodeList.count; ++i) {
          const target = ParticleValues.get(nodeList.targets, i);
          keptOf.set(target, (keptOf.get(target) ?? 0) + 1);
        }
        for (let i = 0; i < nodeList.count; ++i) {
          const target = ParticleValues.get(nodeList.targets, i);
          expect(ParticleValues.get(nodeList.scales!, i)).toBeCloseTo(
            Math.cbrt(byTarget.get(target)! / keptOf.get(target)!),
            5,
          );
        }
        const targets = new Set<number>();
        for (let i = 0; i < nodeList.count; ++i) {
          // a copy of a particle of its tiles: its key that particle's index
          const k = nodeList.keys![i];
          expect(ParticleTiles.getLevelCell(4, tileOf[k], level)).toBe(cell);
          expect(words(nodeList, i)).toEqual(words(whole, k));
          expect(ParticleValues.get(nodeList.targets, i)).toBe(ParticleValues.get(whole.targets, k));
          targets.add(ParticleValues.get(nodeList.targets, i));
        }
        expect(Array.from(node.targets)).toEqual([...targets].sort((a, b) => a - b));
        for (let s = node.base; s < node.base + node.count; ++s) {
          expect(slots.has(s)).toBe(false);
          slots.add(s);
        }
      }
    }
    // the slots: the particles', then the nodes'
    expect(set.count).toBe(base.count + slots.size);
    // its levels left out: its tiles alone
    const tiles = await createParticleTileSetFromPack(stream, { levels: false });
    expect(tiles.levels).toBeUndefined();
    expect([...tiles.tiles.keys()]).toEqual([...base.tiles.keys()]);
    // a pack without: none
    const plain = await openParticlePackStream(ParticlePackSource.fromBytes(writeParticlePack(list, { tiles: 4 })));
    expect(plain.document.tiles!.levels).toBeUndefined();
    expect((await createParticleTileSetFromPack(plain)).levels).toBeUndefined();
  });

  it('reads a blob by slices, a pack with fibers without them; refuses a pack by target as a tile set', async () => {
    const list = createList();
    const tiled = writeParticlePack(list, { tiles: 2 });
    const blob = new Blob([tiled]);
    const set = await createParticleTileSetFromPack(await openParticlePackStream(ParticlePackSource.fromBlob(blob)));
    const tile = await resolve(set.provider.getTile(set.tiles.keys().next().value!, 0));
    expect(tile.count).toBe(set.tiles.values().next().value!.count);

    const byTarget = await openParticlePackStream(ParticlePackSource.fromBytes(writeParticlePack(list)));
    expect(isParticlePackStreamable(byTarget.document)).toBe(false);
    await expect(createParticleTileSetFromPack(byTarget)).rejects.toThrow(/not read by range/);
    // fibers not read (a tile's list none, as a whole list's tiles)
    const fibers = { count: 1, offsets: Int32Array.from([0, 2]), indices: Int32Array.from([1, 5]) };
    const withFibers = writeParticlePack({ ...list, fibers }, { tiles: 2 });
    const fiberStream = await openParticlePackStream(ParticlePackSource.fromBytes(withFibers));
    expect(isParticlePackStreamable(fiberStream.document)).toBe(true);
    expect((await createParticleTileSetFromPack(fiberStream)).count).toBe(list.count);
    await expect(openParticlePackStream(ParticlePackSource.fromBytes(new Uint8Array(100)))).rejects.toThrow(
      /not a particle pack/,
    );
  });
});

describe('ParticlePackSource.batched', () => {
  it('reads ranges asked together near each other as one, at most so many at a time, each its own bytes', async () => {
    const data = Uint8Array.from({ length: 4096 }, (_, i) => i % 251);
    const inner = ParticlePackSource.fromBytes(data);
    let inFlight = 0;
    let most = 0;
    const read = vi.fn(async (offset: number, length: number) => {
      most = Math.max(most, ++inFlight);
      await new Promise((r) => setTimeout(r, 1));
      inFlight -= 1;
      return inner.read(offset, length);
    });
    const source = ParticlePackSource.batched(
      { byteLength: data.byteLength, read },
      { concurrency: 2, gap: 16, maxBytes: 1024 },
    );
    // three side by side (one read), one far (another), two more far apart (two more, at most two in flight)
    const ranges: [number, number][] = [
      [0, 8],
      [8, 8],
      [20, 4],
      [2000, 10],
      [3000, 5],
      [3500, 5],
    ];
    const parts = await Promise.all(ranges.map(([o, l]) => source.read(o, l)));
    for (let k = 0; k < ranges.length; ++k) {
      const [o, l] = ranges[k];
      expect(Array.from(parts[k])).toEqual(Array.from(data.subarray(o, o + l)));
    }
    expect(read).toHaveBeenCalledTimes(4);
    expect(read).toHaveBeenCalledWith(0, 24, expect.anything());
    expect(most).toBeLessThanOrEqual(2);
  });

  it('reads in the order asked (the loader its priority), not by offset; a range aborted pending not read', async () => {
    const data = Uint8Array.from({ length: 4096 }, (_, i) => i % 251);
    const inner = ParticlePackSource.fromBytes(data);
    const order: number[] = [];
    const read = vi.fn(async (offset: number, length: number) => {
      order.push(offset);
      await new Promise((r) => setTimeout(r, 1));
      return inner.read(offset, length);
    });
    const source = ParticlePackSource.batched({ byteLength: data.byteLength, read }, { concurrency: 1, gap: 16 });
    const dropped = new AbortController();
    const reads = [
      source.read(3000, 8),
      source.read(1000, 8, dropped.signal),
      source.read(10, 8),
      // near the first: merged into its read
      source.read(3010, 8),
    ];
    dropped.abort(new Error('no longer wanted'));
    const [far, aborted, near, merged] = await Promise.allSettled(reads);
    expect(far.status).toBe('fulfilled');
    expect(near.status).toBe('fulfilled');
    expect(merged.status).toBe('fulfilled');
    expect(aborted.status).toBe('rejected');
    expect(String((aborted as PromiseRejectedResult).reason)).toContain('no longer wanted');
    expect(order).toEqual([3000, 10]);
  });

  it('aborts a read in flight where every range of it is aborted, not while one is still wanted', async () => {
    const data = new Uint8Array(4096);
    const signals: AbortSignal[] = [];
    const read = vi.fn(
      (_offset: number, length: number, signal?: AbortSignal) =>
        new Promise<Uint8Array>((resolve, reject) => {
          signals.push(signal!);
          const t = setTimeout(() => resolve(new Uint8Array(length)), 20);
          signal?.addEventListener('abort', () => {
            clearTimeout(t);
            reject(signal.reason);
          });
        }),
    );
    const source = ParticlePackSource.batched({ byteLength: data.byteLength, read }, { gap: 16 });
    const a = new AbortController();
    const b = new AbortController();
    const reads = [source.read(0, 8, a.signal), source.read(8, 8, b.signal)];
    await new Promise((r) => setTimeout(r, 1));
    expect(read).toHaveBeenCalledTimes(1);
    a.abort(new Error('a let go'));
    expect(signals[0].aborted).toBe(false);
    b.abort(new Error('b let go'));
    expect(signals[0].aborted).toBe(true);
    const settled = await Promise.allSettled(reads);
    expect(settled.map((s) => s.status)).toEqual(['rejected', 'rejected']);
  });
});

describe('a nested particle pack read by range', () => {
  /** a template of 40,000 particles (two tiles an axis), placed by two copies */
  function createNestedList(): ParticleList {
    let seed = 5;
    const random = () => ((seed = (seed * 1103515245 + 12345) % 2 ** 31) / 2 ** 31) * 200 - 100;
    const n = 40_000;
    const template: ParticleList = {
      label: 'spread',
      count: n,
      coordinates: Float32Array.from({ length: n * 3 }, random),
      rotations: Float32Array.from({ length: n * 4 }, (_, k) => (k % 4 === 3 ? 1 : 0)),
      targets: Int32Array.from({ length: n }, (_, i) => i % 3),
      targetInfo: new Map([
        [0, { name: 'a' }],
        [1, { name: 'b' }],
        [2, { name: 'c' }],
      ]),
      getParticleLabel: (i) => `${i}`,
      sourceData: { kind: 'test', name: 'test', data: undefined },
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    return {
      label: 'nested',
      count: 2,
      coordinates: Float32Array.of(0, 0, 0, 300, 0, 0),
      targets: new Int32Array(2),
      targetInfo: new Map([[0, { name: 'copy' }]]),
      targetMapping: new Map<number, ParticleTarget>([[0, { kind: 'particles', particles: template }]]),
      getParticleLabel: (i) => `${i}`,
      sourceData: { kind: 'test', name: 'test', data: undefined },
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
  }

  it("reads its templates' runs, directories and levels, their particles by range when asked", async () => {
    const bytes = writeParticlePack(createNestedList());
    const whole = createParticleListFromPack(readParticlePack(bytes));
    const stream = await openParticlePackStream(ParticlePackSource.fromBytes(bytes));
    expect(isParticlePackListStreamable(stream.document)).toBe(true);
    const list = await openParticleListFromPackStream(stream);
    const [tw, ts] = [whole, list].map((l) => ParticleNesting.getTemplate(l, 0)!.particles);
    // the template deferred: its arrays empty (their kinds), its counts and runs as read whole
    const deferred = Particle.getDeferred(ts)!;
    expect(deferred).toBeDefined();
    expect(Particle.getDeferred(tw)).toBeUndefined();
    expect(ts.count).toBe(tw.count);
    expect(GridPositions.is(ts.coordinates) && GridPositions.count(ts.coordinates)).toBe(0);
    expect(ts.targets).toEqual(tw.targets);
    expect(ParticleTiles.getDirectory(ts)).toBeDefined();
    // its levels stored, as the whole read's
    const [lw, ls] = [whole, list].map((l, i) => ParticleNesting.getTemplateLevels(l, [tw, ts][i]));
    expect(ls.map((l) => Array.from(l.indices))).toEqual(lw.map((l) => Array.from(l.indices)));
    // a range of its particles read: as the whole read's
    const first = 8192;
    const count = 8192;
    const read = await deferred.read(first, count);
    expect(Array.from(read.positions)).toEqual(
      Array.from((tw.coordinates as GridPositions).data.subarray(first * 2, (first + count) * 2)),
    );
    expect(Array.from(read.rotations!)).toEqual(Array.from(tw.rotations!.subarray(first * 4, (first + count) * 4)));
    // what reads its particles refuses it (a subset gathered)
    expect(() => Particle.subset(ts, Int32Array.of(0, 1))).toThrow(/read by range/);
    expect(() => ParticleNesting.flatten(list)).toThrow(/read by range/);
    expect(() => writeParticlePack(list)).toThrow(/read by range/);
    expect(Vec3.distance(Particle.getBoundary(ts).sphere.center, Particle.getBoundary(tw).sphere.center)).toBe(0);
  });

  it("bounds a loci of its copies held in part by its template's tiles, not its particles", async () => {
    const bytes = writeParticlePack(createNestedList());
    const whole = createParticleListFromPack(readParticlePack(bytes));
    const list = await openParticleListFromPackStream(
      await openParticlePackStream(ParticlePackSource.fromBytes(bytes)),
    );
    const [tw, ts] = [whole, list].map((l) => ParticleNesting.getTemplate(l, 0)!.particles);
    /** the spheres of `loci` of the whole read and the read by range: the latter holds the former */
    const spheres = (of: (l: ParticleList, t: ParticleList) => Particle.Loci) => {
      const [w, s] = [whole, list].map((l, i) => Particle.getBoundingSphere(of(l, [tw, ts][i]), Sphere3D()));
      expect(Vec3.distance(w.center, s.center) + w.radius).toBeLessThanOrEqual(s.radius + 1e-3);
      return [w, s];
    };
    // a page's particles in both copies, a particle of one copy (a set too small to bound by its own sphere)
    const [w, s] = spheres((l, t) => ParticleNesting.selectInTemplate(l, t, (j) => j < 8192));
    expect(s.radius).toBeLessThan(w.radius * 1.5);
    spheres((l) => ParticleNesting.ofInstance(l, 1, 5));
    // a copy whole: its template's boundary, as read whole
    const [cw, cs] = spheres((l) => ParticleNesting.ofCopy(l, 0));
    expect(cs.radius).toBeCloseTo(cw.radius, 3);
  });

  it("keeps its templates' pages read in one cache: a page read once, a loci bounded by its particles where kept", async () => {
    const bytes = writeParticlePack(createNestedList());
    const whole = createParticleListFromPack(readParticlePack(bytes));
    const tw = ParticleNesting.getTemplate(whole, 0)!.particles;
    /** the list read by range from a source counting its reads, its cache of `pageCacheBytes` */
    const open = async (pageCacheBytes?: number) => {
      const inner = ParticlePackSource.fromBytes(bytes);
      const counter = { reads: 0 };
      const source = {
        byteLength: inner.byteLength,
        read: (offset: number, length: number) => {
          counter.reads += 1;
          return inner.read(offset, length);
        },
      };
      const list = await openParticleListFromPackStream(await openParticlePackStream(source), { pageCacheBytes });
      counter.reads = 0;
      return { list, counter, deferred: Particle.getDeferred(ParticleNesting.getTemplate(list, 0)!.particles)! };
    };
    const { list, counter, deferred } = await open();
    // a particle of the second copy: by its tile's box, no read begun
    const loci = ParticleNesting.ofInstance(list, 1, 5);
    const exact = Particle.getBoundingSphere(ParticleNesting.ofInstance(whole, 1, 5), Sphere3D());
    const byTile = Particle.getBoundingSphere(loci, Sphere3D());
    expect(counter.reads).toBe(0);
    expect(deferred.getPage(0)).toBeUndefined();
    // a particle of a stored level (a coarse copy's click): by its position there, its page not read
    const stored = Particle.getStoredLevels(ParticleNesting.getTemplate(list, 0)!.particles)!;
    const coarse = stored[0].keys.find((k) => k >> 13 === 1)!;
    const inLevel = Particle.getBoundingSphere(ParticleNesting.ofInstance(list, 1, coarse), Sphere3D());
    const inWhole = Particle.getBoundingSphere(ParticleNesting.ofInstance(whole, 1, coarse), Sphere3D());
    expect(Vec3.distance(inLevel.center, inWhole.center)).toBeLessThan(1e-3);
    expect(inLevel.radius).toBeCloseTo(inWhole.radius, 3);
    expect(counter.reads).toBe(0);
    expect(deferred.getPage(1)).toBeUndefined();
    // its page read (its positions and rotations: a read each), once while in flight, kept: by its position then
    const [p0, again] = await Promise.all([deferred.readPage(0), deferred.readPage(0)]);
    expect(again).toBe(p0);
    expect(counter.reads).toBe(2);
    expect(Array.from(p0.positions)).toEqual(Array.from((tw.coordinates as GridPositions).data.subarray(0, 8192 * 2)));
    expect(deferred.getPage(0)).toBe(p0);
    await deferred.readPage(0);
    expect(counter.reads).toBe(2);
    const byPosition = Particle.getBoundingSphere(loci, Sphere3D());
    expect(Vec3.distance(byPosition.center, exact.center)).toBeLessThan(1e-3);
    expect(byPosition.radius).toBeCloseTo(exact.radius, 3);
    expect(byPosition.radius).toBeLessThan(byTile.radius / 10);

    // a cache of a page's bytes: the next page read lets the first go, its loci by its tile again
    const small = await open(p0.positions.byteLength + (p0.rotations?.byteLength ?? 0));
    await small.deferred.readPage(0);
    await small.deferred.readPage(1);
    expect(small.deferred.getPage(0)).toBeUndefined();
    expect(small.deferred.getPage(1)).toBeDefined();
    const later = Particle.getBoundingSphere(ParticleNesting.ofInstance(small.list, 1, 5), Sphere3D());
    expect(later.radius).toBeCloseTo(byTile.radius, 3);
  });
});
