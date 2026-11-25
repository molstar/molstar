/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { Cell } from '@molstar/core/math/geometry/spacegroup/cell';
import { RuntimeContext } from '@molstar/core/task';
import { CIF } from '@molstar/io/reader/cif';
import { CustomProperties } from '../../../model/custom-property.js';
import {
  Particle,
  type ParticleList,
  ParticleRotations,
  ParticleSites,
  ParticleTarget,
  ParticleValues,
} from '../../../model/particles/particle-list.js';
import { ParticleTileCatalog, ParticleTileGrid, ParticleTiles } from '../../../model/particles/particle-tiles.js';
import { mortonOrder3d } from '@molstar/core/data/util/hash-functions';
import { ParticleTileView } from '../../../model/particles/particle-tile-view.js';
import { ParticleSlots } from '../../../model/particles/particle-slots.js';
import { ParticleNesting } from '../../../model/particles/particle-nesting.js';
import { ParticleDensityGrid, ParticleTileOpaque } from '../../../model/particles/particle-tile-opaque.js';
import { createParticleListFromMmcifAssembly } from '../mmcif.js';
import {
  createParticleListFromPack,
  isParticlePack,
  ParticlePackFormat,
  parseParticlePackDocument,
  readParticlePack,
  writeParticlePack,
} from '../pack.js';

/** `actual` within `eps` of `expected`, value by value (quantized sites: within half the quantization step) */
function expectClose(actual: ArrayLike<number>, expected: ArrayLike<number>, eps: number) {
  expect(actual.length).toBe(expected.length);
  for (let k = 0; k < expected.length; ++k) expect(Math.abs(actual[k] - expected[k])).toBeLessThanOrEqual(eps);
}

/**
 * Seven particles of three targets, interleaved (the pack sorts them by target): rotations, a radius a target (or
 * differing in rounding),
 * entities and compartments by target, two targets with sites (one with chains), a fiber, an attribute, a cell
 */
function createList(options: { roundedRadii?: boolean } = {}): ParticleList {
  const targets = Int32Array.from([5, 2, 5, 9, 2, 5, 9]);
  const n = targets.length;
  const coordinates = new Float32Array(n * 3).map((_, k) => k * 0.5 - 3);
  const rotations = new Float32Array(n * 4);
  for (let i = 0; i < n; ++i) {
    const a = i * 0.3;
    rotations.set([Math.sin(a), 0, 0, Math.cos(a)], i * 4);
  }
  const byTarget = new Map([
    [5, { entity: 0, compartment: 1, radius: 4 }],
    [2, { entity: 1, compartment: -1, radius: 2 }],
    [9, { entity: 1, compartment: 0, radius: 7 }],
  ]);
  const entities = Int32Array.from(targets, (t) => byTarget.get(t)!.entity);
  const compartments = Int32Array.from(targets, (t) => byTarget.get(t)!.compartment);
  const radii = Float32Array.from(targets, (t, i) => byTarget.get(t)!.radius + (options.roundedRadii ? i * 1e-6 : 0));
  const sites5 = {
    count: 3,
    positions: Float32Array.from([0, 0, 0, 1, 0, 0, 0, 2, 0]),
    radii: Float32Array.from([1.7, 1.55, 1.52]),
    label: 'A',
    chains: { offsets: Int32Array.from([0, 2, 3]), labels: ['A', 'B'] },
  };
  const sites9 = { count: 1, positions: Float32Array.from([3, 3, 3]), radii: Float32Array.from([2]) };
  const score = Float32Array.from([0.5, 0.25, 1, 0, 0.75, 0.125, 0.875]);
  return {
    label: 'pack test',
    entryId: 'TEST',
    count: n,
    cell: Cell.create(Vec3.create(100, 200, 300), Vec3.create(Math.PI / 2, Math.PI / 2, Math.PI / 3)),
    keys: Int32Array.from({ length: n }, (_, i) => 100 + i),
    targets,
    targetInfo: new Map([
      [5, { name: 'five', entity: 0 }],
      [2, { name: 'two', entity: 1 }],
      [9, { name: 'nine', entity: 1 }],
    ]),
    targetMapping: new Map<number, ParticleTarget>([
      [5, { kind: 'sites', sites: sites5, center: Vec3.create(0.25, 0.5, 0) }],
      [9, { kind: 'sites', sites: sites9 }],
    ]),
    entities,
    entityInfo: new Map([
      [0, { name: 'ent0', function: 'Lipid' }],
      [1, { name: 'ent1' }],
    ]),
    compartments,
    compartmentInfo: new Map([
      [0, { name: 'root.a' }],
      [1, { name: 'root.b' }],
    ]),
    coordinates,
    rotations,
    radii,
    // target 5's particles 0, 2, 5 as a fiber, in that order
    fibers: { count: 1, offsets: Int32Array.from([0, 3]), indices: Int32Array.from([0, 2, 5]) },
    attributes: new Map([['score', { column: Column.ofFloatArray(score), label: 'Score', min: 0, max: 1 }]]),
    getParticleLabel: (i) => `p${i}`,
    sourceData: { kind: 'test', name: 'test', data: undefined },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

/**
 * 300 particles of three targets spread over a box (an entity and a radius a target), a fiber, their indices an
 * attribute (`source`: a particle's in the list, through any order)
 */
function createSpread(): ParticleList {
  const n = 300;
  const coordinates = new Float32Array(n * 3);
  for (let i = 0; i < n; ++i) coordinates.set([(i * 37) % 101, (i * 53) % 97, (i * 11) % 89], i * 3);
  const targets = Int32Array.from({ length: n }, (_, i) => [3, 1, 2][i % 3]);
  return {
    label: 'spread',
    count: n,
    targets,
    targetInfo: new Map([
      [1, { name: 'one', entity: 0 }],
      [2, { name: 'two', entity: 1 }],
      [3, { name: 'three', entity: 1 }],
    ]),
    entities: Int32Array.from(targets, (t) => (t === 1 ? 0 : 1)),
    entityInfo: new Map([
      [0, { name: 'e0' }],
      [1, { name: 'e1' }],
    ]),
    coordinates,
    rotations: new Float32Array(n * 4).map((_, k) => (k % 4 === 3 ? 1 : 0)),
    radii: Float32Array.from(targets, (t) => t * 2),
    fibers: { count: 1, offsets: Int32Array.from([0, 4]), indices: Int32Array.from([7, 100, 3, 250]) },
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

/** `createSpread`'s list, its targets of sites (a volume each: the tiles' opacities') */
function createSpreadWithSites(): ParticleList {
  const list = createSpread();
  const sites = (n: number, r: number): ParticleTarget => ({
    kind: 'sites',
    sites: {
      count: n,
      positions: Float32Array.from({ length: n * 3 }, (_, k) => k % 3),
      radii: new Float32Array(n).fill(r),
    },
  });
  return {
    ...list,
    targetMapping: new Map([
      [1, sites(4, 3)],
      [2, sites(2, 5)],
      [3, sites(1, 8)],
    ]),
  };
}

/** The particles' order in a pack of `list`: by target in the order of their first particle, each target's in order */
function getPackOrder(list: ParticleList) {
  const order: number[] = [];
  const seen: number[] = [];
  for (let i = 0; i < list.count; ++i) if (!seen.includes(list.targets[i])) seen.push(list.targets[i]);
  for (const t of seen) for (let i = 0; i < list.count; ++i) if (list.targets[i] === t) order.push(i);
  return order;
}

/** Particle `p` of `read` (a pack's: on a grid) at particle `i` of `list`, within half the grid's step */
function expectPosition(read: ParticleList, p: number, list: ParticleList, i: number) {
  if (!GridPositions.is(read.coordinates)) throw new Error('positions on a grid expected');
  const a = GridPositions.get(read.coordinates, p, Vec3());
  const b = GridPositions.get(list.coordinates, i, Vec3());
  for (let k = 0; k < 3; ++k) expect(Math.abs(a[k] - b[k])).toBeLessThanOrEqual(read.coordinates.step / 2 + 1e-6);
}

describe('particle pack', () => {
  it('copies a pack at an unaligned offset only where small, else refuses (no copy of its size again)', () => {
    const bytes = writeParticlePack(createList());
    const shifted = new Uint8Array(bytes.byteLength + 1);
    shifted.set(bytes, 1);
    const unaligned = shifted.subarray(1);
    expect(createParticleListFromPack(readParticlePack(unaligned)).count).toBe(
      createParticleListFromPack(readParticlePack(bytes)).count,
    );
    expect(() => readParticlePack(unaligned, { maxCopyBytes: bytes.byteLength - 1 })).toThrow(/not copied/);
    // aligned: never copied, whatever its size
    expect(readParticlePack(bytes, { maxCopyBytes: 0 }).buffer).toBe(bytes.buffer);
  });

  it('reads no block whose sizes disagree, naming it: its bytes and items, a particle block and the count', () => {
    const bytes = writeParticlePack(createList());
    const { document } = readParticlePack(bytes);
    const parse = (change: (blocks: { name: string; count: number; byteLength: number }[]) => void) => {
      const doc = JSON.parse(JSON.stringify(document));
      change(doc.blocks);
      return () => parseParticlePackDocument(new TextEncoder().encode(JSON.stringify(doc)), bytes.byteLength);
    };
    expect(parse(() => {})).not.toThrow();
    // its items more than its bytes hold: the next block's read
    expect(
      parse((blocks) => {
        blocks.find((b) => b.name === 'rotation')!.count += 1;
      }),
    ).toThrow(/block 'rotation' of \d+ bytes, not its/);
    // its bytes and items agree, of another count than the particles'
    expect(
      parse((blocks) => {
        const b = blocks.find((b) => b.name === 'rotation')!;
        b.byteLength -= b.byteLength / b.count;
        b.count -= 1;
      }),
    ).toThrow(/block 'rotation' of \d+ items, not the pack's/);
  });

  it('reads back the list it was written from, its particles sorted by target', () => {
    const list = createList();
    const bytes = writeParticlePack(list);
    expect(isParticlePack(bytes)).toBe(true);
    const read = createParticleListFromPack(readParticlePack(bytes));
    const order = getPackOrder(list);

    expect(read.count).toBe(list.count);
    // no keys: a particle's key is its index in the pack
    expect(read.keys).toBeUndefined();
    expect(read.label).toBe('pack test');
    expect(read.entryId).toBe('TEST');
    expect(Array.from(read.cell!.size)).toEqual([100, 200, 300]);
    for (let k = 0; k < 3; ++k) expect(read.cell!.anglesInRadians[k]).toBeCloseTo(list.cell!.anglesInRadians[k], 12);
    for (let p = 0; p < read.count; ++p) {
      const i = order[p];
      expect(ParticleValues.get(read.targets!, p)).toBe(list.targets[i]);
      expect(ParticleValues.get(read.entities!, p)).toBe(list.entities![i]);
      expect(ParticleValues.get(read.compartments!, p)).toBe(list.compartments![i]);
      expect(ParticleValues.get(read.radii!, p)).toBe(list.radii![i]);
      expectPosition(read, p, list, i);
      // snorm16 (`ParticleRotations`)
      expect(Array.from(read.rotations!.subarray(p * 4, p * 4 + 4))).toEqual(
        Array.from(ParticleRotations.toSnorm16(list.rotations!).subarray(i * 4, i * 4 + 4)),
      );
      expect(read.attributes!.get('score')!.column.value(p)).toBe(list.attributes!.get('score')!.column.value(i));
    }
    expect(read.attributes!.get('score')).toMatchObject({ label: 'Score', min: 0, max: 1 });
    // a target's info and its compartment, stored once (its particles' all)
    expect([...read.targetInfo]).toEqual(
      [...list.targetInfo].map(([t, info]) => [t, { ...info, compartment: { 5: 1, 2: -1, 9: 0 }[t] }]),
    );
    expect([...read.entityInfo!]).toEqual([...list.entityInfo!]);
    expect([...read.compartmentInfo!]).toEqual([...list.compartmentInfo!]);
    // the fiber's particles at their new positions, in their order
    expect(read.fibers!.count).toBe(1);
    expect(Array.from(read.fibers!.indices)).toEqual([0, 1, 2]);
    expect(Array.from(read.fibers!.indices, (p) => order[p])).toEqual([0, 2, 5]);

    // the sites, with their chains and centers (given, or the centroid)
    const five = read.targetMapping!.get(5)!;
    const nine = read.targetMapping!.get(9)!;
    expect(read.targetMapping!.has(2)).toBe(false);
    if (five.kind !== 'sites' || nine.kind !== 'sites') throw new Error('sites expected');
    // quantized: within half a step (the extent 2 over 32767)
    expectClose(five.sites.positions, [0, 0, 0, 1, 0, 0, 0, 2, 0], 1 / 32767);
    expect(Array.from(five.sites.radii)).toEqual(Array.from(Float32Array.from([1.7, 1.55, 1.52])));
    expect(five.sites.label).toBe('A');
    expect(Array.from(five.sites.chains!.offsets)).toEqual([0, 2, 3]);
    expect(five.sites.chains!.labels).toEqual(['A', 'B']);
    expect(Array.from(five.center!)).toEqual([0.25, 0.5, 0]);
    expect(Array.from(nine.sites.positions)).toEqual([3, 3, 3]);
    expect(Array.from(ParticleTarget.getCenter(nine)!)).toEqual([3, 3, 3]);
    expect(nine.sites.chains).toBeUndefined();

    expect(read.getParticleLabel(0)).toBe('Particle 1 · ent0 · target five');
    expect(ParticlePackFormat.is(read.sourceData)).toBe(true);
  });

  it('quantizes the sites a target at a time: within half a step of its extent, radii by a table of types', () => {
    const list = createList();
    // a fiber's extent: 3000 Å, a step of 3000 / 65534 Å
    const count = 1000;
    const positions = new Float32Array(count * 3).map((_, k) => ((k * 7919) % 3001) - 1500 + (k % 3) * 0.37);
    const radii = Float32Array.from({ length: count }, (_, k) => [1.7, 1.55, 1.52, 1.8][k % 4]);
    const sites = { count, positions, radii };
    const pack = readParticlePack(
      writeParticlePack({
        ...list,
        targetMapping: new Map([...list.targetMapping!, [2, { kind: 'sites', sites } as ParticleTarget]]),
      }),
    );
    // 8 bytes a site, the radii's types a table of their distinct values
    const block = pack.document.blocks.find((b) => b.name === 'site')!;
    expect(block.byteLength).toBe((3 + 1 + count) * 8);
    expect([...pack.document.siteRadii].sort()).toEqual([1.52, 1.55, 1.7, 1.8, 2].map(Math.fround).sort());
    const read = createParticleListFromPack(pack).targetMapping!.get(2)!;
    if (read.kind !== 'sites') throw new Error('sites expected');
    const step = read.sites.quantized!.scale;
    expect(step).toBeCloseTo(3000 / 2 / 32767, 3);
    expectClose(read.sites.positions, positions, step / 2 + 1e-4);
    expect(Array.from(read.sites.radii)).toEqual(Array.from(radii));
    // a site at a time as the arrays decoded
    const p = Vec3();
    for (const i of [0, 17, count - 1]) {
      expect(Array.from(ParticleSites.getPosition(read.sites, i, p), Math.fround)).toEqual(
        Array.from(read.sites.positions.subarray(i * 3, i * 3 + 3)),
      );
      expect(ParticleSites.getRadius(read.sites, i)).toBe(radii[i]);
    }
  });

  it("keeps the sites' elements by type, a type a distinct radius and element (S and P share a radius)", () => {
    const list = createList();
    // C, N, O, S, P: S and P of one radius
    const count = 10;
    const positions = new Float32Array(count * 3).map((_, k) => k * 0.5);
    const elements = Uint8Array.from({ length: count }, (_, k) => [6, 7, 8, 16, 15][k % 5]);
    const radii = Float32Array.from(elements, (e) => ({ 6: 1.7, 7: 1.55, 8: 1.52, 16: 1.8, 15: 1.8 })[e]!);
    const withElements = (t: ParticleTarget): ParticleTarget =>
      t.kind === 'sites' ? { ...t, sites: { ...t.sites, elements: new Uint8Array(t.sites.count).fill(6) } } : t;
    const targetMapping = new Map<number, ParticleTarget>([
      ...[...list.targetMapping!].map(([id, t]) => [id, withElements(t)] as const),
      [2, { kind: 'sites', sites: { count, positions, radii, elements } }],
    ]);
    const pack = readParticlePack(writeParticlePack({ ...list, targetMapping }));
    const types = pack.document.siteRadii.map((r, k) => [r, pack.document.siteElements![k]]);
    // the other targets' radii 1.55 and 1.52 with carbon, 2 with carbon: types of their own
    expect(
      types
        .filter(([r]) => r === Math.fround(1.8))
        .map(([, e]) => e)
        .sort(),
    ).toEqual([15, 16]);
    const read = createParticleListFromPack(pack).targetMapping!.get(2)!;
    if (read.kind !== 'sites') throw new Error('sites expected');
    expect(ParticleSites.hasElements(read.sites)).toBe(true);
    expect(Array.from(read.sites.elements!)).toEqual(Array.from(elements));
    expect(Array.from(read.sites.radii)).toEqual(Array.from(radii));
    for (let i = 0; i < count; ++i) expect(ParticleSites.getElement(read.sites, i)).toBe(elements[i]);
  });

  it("writes no elements where a target's sites have none, and reads a version 4 pack without them", () => {
    const bytes = writeParticlePack(createList());
    expect(readParticlePack(bytes).document.siteElements).toBeUndefined();
    // the version before elements: its document the same but for them
    new DataView(bytes.buffer).setUint32(4, 4, true);
    const read = createParticleListFromPack(readParticlePack(bytes)).targetMapping!.get(5)!;
    if (read.kind !== 'sites') throw new Error('sites expected');
    expect(ParticleSites.hasElements(read.sites)).toBe(false);
    expect(read.sites.elements).toBeUndefined();
    expect(ParticleSites.getElement(read.sites, 0)).toBe(0);
  });

  it("keeps each chain's extremes of its quantized positions (one for a target without chains)", () => {
    const pack = readParticlePack(writeParticlePack(createList()));
    const read = createParticleListFromPack(pack);
    for (const id of [5, 9]) {
      const target = read.targetMapping!.get(id)!;
      if (target.kind !== 'sites') throw new Error('sites expected');
      const q = target.sites.quantized!;
      const runs = target.sites.chains?.offsets ?? Int32Array.of(0, target.sites.count);
      expect(q.extremes!.length).toBe((runs.length - 1) * 6);
      for (let c = 0; c + 1 < runs.length; ++c) {
        const expected = ParticleSites.getQuantizedExtremes(q.data, runs[c], runs[c + 1], new Int16Array(6));
        expect(Array.from(q.extremes!.subarray(c * 6, c * 6 + 6))).toEqual(Array.from(expected));
      }
    }
  });

  it('reads its arrays as views of the bytes, each block aligned to 256 bytes', () => {
    const bytes = writeParticlePack(createList());
    const pack = readParticlePack(bytes);
    const read = createParticleListFromPack(pack);
    for (const b of pack.document.blocks) expect(b.byteOffset % 256).toBe(0);
    expect(GridPositions.is(read.coordinates) && read.coordinates.data.buffer).toBe(bytes.buffer);
    expect(read.rotations!.buffer).toBe(bytes.buffer);
    const five = read.targetMapping!.get(5)!;
    if (five.kind !== 'sites') throw new Error('sites expected');
    expect(five.sites.quantized!.data.buffer).toBe(bytes.buffer);
    expect(read.fibers!.indices.buffer).toBe(bytes.buffer);
  });

  it("keeps a target's radius once, the largest where its particles' differ in rounding", () => {
    const list = createList({ roundedRadii: true });
    const pack = readParticlePack(writeParticlePack(list));
    expect(pack.document.targets.map((t) => t.radius)).toEqual([4 + 5e-6, 2 + 4e-6, 7 + 6e-6].map(Math.fround));
    const read = createParticleListFromPack(pack);
    for (const t of pack.document.targets) {
      for (let p = t.particleOffset; p < t.particleOffset + t.particleCount; ++p)
        expect(ParticleValues.get(read.radii!, p)).toBe(t.radius);
    }
    // what a list has per particle and the pack per target: not stored per particle
    const names = pack.document.blocks.map((b) => b.name);
    for (const name of ['key', 'target', 'entity', 'compartment', 'radius']) expect(names).not.toContain(name);
  });

  it("keeps a target's bounds, rounded outward: read without a pass, computed alike where a pack has none", () => {
    const list = createList();
    const pack = readParticlePack(writeParticlePack(list));
    const read = createParticleListFromPack(pack);
    const position = Vec3();
    for (const t of pack.document.targets) {
      const { min, max, distance } = t.bounds!;
      const start = t.particleOffset;
      const end = start + t.particleCount;
      const b = Particle.getRangeBounds(read, start, end)!;
      expect(Array.from(b.min)).toEqual(min);
      expect(Array.from(b.max)).toEqual(max);
      expect(b.distance).toBe(distance);
      // what it holds: every particle in the box, within the distance of its center
      const center = Vec3.scale(Vec3(), Vec3.add(Vec3(), b.min, b.max), 0.5);
      for (let p = start; p < end; ++p) {
        GridPositions.get(read.coordinates, p, position);
        for (let k = 0; k < 3; ++k) {
          expect(position[k]).toBeGreaterThanOrEqual(min[k]);
          expect(position[k]).toBeLessThanOrEqual(max[k]);
        }
        expect(Vec3.distance(position, center)).toBeLessThanOrEqual(distance);
      }
      // a pack without them: computed from the positions, exact (within the rounding of the stored ones)
      const doc = { ...pack.document, targets: pack.document.targets.map(({ bounds: _, ...rest }) => rest) };
      const computed = Particle.getRangeBounds(createParticleListFromPack({ ...pack, document: doc }), start, end)!;
      expectClose(computed.min, min, 2e-3);
      expectClose(computed.max, max, 2e-3);
      expect(computed.distance).toBeLessThanOrEqual(distance);
      expect(computed.distance).toBeGreaterThan(distance - 3e-3);
    }
    // not a target's range: none; a moved list's computed again, by its positions
    const first = pack.document.targets[0];
    expect(
      Particle.getRangeBounds(read, first.particleOffset, first.particleOffset + first.particleCount - 1),
    ).toBeUndefined();
    const moved = Particle.transform(read, Mat4.fromTranslation(Mat4(), Vec3.create(10, 0, 0)));
    const end = first.particleOffset + first.particleCount;
    const b = Particle.getRangeBounds(moved, first.particleOffset, end)!;
    expect(b.min[0]).toBeCloseTo(Particle.getRangeBounds(read, first.particleOffset, end)!.min[0] + 10, 2);
  });

  it('refuses a target whose particles differ in what it stores once (an entity, a compartment)', () => {
    const list = createList();
    const entities = Int32Array.from(list.entities!);
    entities[2] = 1;
    expect(() => writeParticlePack({ ...list, entities })).toThrow(/target 5's particles differ in entity/);
  });

  it('reads a pack at an offset that does not keep its blocks aligned, copied once', () => {
    const bytes = writeParticlePack(createList());
    const shifted = new Uint8Array(bytes.length + 1);
    shifted.set(bytes, 1);
    const read = createParticleListFromPack(readParticlePack(shifted.subarray(1)));
    expect(read.count).toBe(7);
    expect(GridPositions.get(read.coordinates, 0, Vec3())[0]).toBeCloseTo(-3, 4);
  });

  it('refuses what is not a pack, or of another version', () => {
    expect(() => readParticlePack(new Uint8Array(128))).toThrow(/MSPK/);
    const bytes = writeParticlePack(createList());
    new DataView(bytes.buffer).setUint32(4, 99, true);
    expect(() => readParticlePack(bytes)).toThrow(/version 99/);
    // a version 5 pack (sorted by target, as 6 without tiles) read as it is
    new DataView(bytes.buffer).setUint32(4, 5, true);
    expect(createParticleListFromPack(readParticlePack(bytes)).count).toBe(7);
  });

  it("writes a tile-ordered pack's tiles' opacities, its tile sets' as computed (a coarser one's its cells' blocks)", async () => {
    const list = createSpreadWithSites();
    const read = createParticleListFromPack(readParticlePack(writeParticlePack(list, { tiles: 4 })));
    const plain: ParticleList = {
      ...read,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    expect(ParticleTiles.getDirectory(read)!.opacity!.blocks).toBe(4);
    for (const [d, blocks] of [
      [4, 4],
      [2, 8],
    ]) {
      const fromPack = ParticleTiles.fromList(read, [d, d, d]);
      const computed = await ParticleTileOpaque.computeOpacities(ParticleTiles.fromList(plain, [d, d, d]), {
        blocks,
      }).run();
      let any = false;
      for (const [id, tile] of fromPack.tiles) {
        const other = computed.tiles.get(id)!.opacity;
        expect(tile.opacity!.dimensions).toEqual([blocks, blocks, blocks]);
        if (!other) {
          expect(tile.opacity!.density.every((v) => v === 0)).toBe(true);
          continue;
        }
        any = true;
        if (d === 4) {
          for (let k = 0; k < other.density.length; ++k)
            expect(tile.opacity!.density[k]).toBeCloseTo(other.density[k], 6);
        } else {
          // its cells' blocks: a molecule spread over its own cell's blocks alone, so block by block not the coarse
          // tile's own computation, but the same cross-sections in all
          const sum = (a: Float32Array) => a.reduce((x, v) => x + v, 0);
          expect(sum(tile.opacity!.density)).toBeCloseTo(sum(other.density), 5);
        }
      }
      expect(any).toBe(true);
    }
    // none written where asked: none read
    const none = createParticleListFromPack(readParticlePack(writeParticlePack(list, { tiles: 4, opacityBlocks: 0 })));
    expect(ParticleTiles.getDirectory(none)!.opacity).toBeUndefined();
    expect([...ParticleTiles.fromList(none, [4, 4, 4]).tiles.values()].every((t) => !t.opacity)).toBe(true);
  });

  it("writes a pack by target's density grid, as computed from it read without one; none tile-ordered", () => {
    const list = createSpreadWithSites();
    const read = createParticleListFromPack(readParticlePack(writeParticlePack(list, { densityBlocks: 4 })));
    const grid = ParticleDensityGrid.get(read)!;
    const plain: ParticleList = {
      ...read,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    const computed = ParticleDensityGrid.compute(plain, { longest: 4 })!;
    expect(grid.dimensions).toEqual(computed.dimensions);
    expectClose(grid.min, computed.min, 1e-6);
    expectClose(grid.max, computed.max, 1e-6);
    expectClose(grid.density, computed.density, 1e-9);
    expect(grid.density.some((d) => d > 0)).toBe(true);
    // 64 along the longest side by default
    const byDefault = createParticleListFromPack(readParticlePack(writeParticlePack(list)));
    expect(Math.max(...ParticleDensityGrid.get(byDefault)!.dimensions)).toBe(64);
    // none where asked, none tile-ordered (its tiles' opacities are)
    expect(
      ParticleDensityGrid.get(
        createParticleListFromPack(readParticlePack(writeParticlePack(list, { densityBlocks: 0 }))),
      ),
    ).toBeUndefined();
    const tiled = readParticlePack(writeParticlePack(list, { tiles: 4 }));
    expect(tiled.document.density).toBeUndefined();
    expect(ParticleDensityGrid.get(createParticleListFromPack(tiled))).toBeUndefined();
  });

  it("draws a tile-ordered pack's list whole as a slot list of its tiles: a slot its index, a tile a view of its range", () => {
    const list = createSpreadWithSites();
    const read = createParticleListFromPack(readParticlePack(writeParticlePack(list, { tiles: 4 })));
    const whole = ParticleTileView.ofWhole(read)!;
    expect(ParticleTileView.ofWhole(read)).toBe(whole);
    expect(whole.count).toBe(read.count);
    expect(ParticleSlots.getLiveCount(whole)).toBe(read.count);
    const p = Vec3();
    const q = Vec3();
    for (const page of whole.pages!.pages) {
      expect((page.list.coordinates as GridPositions).data.buffer).toBe(
        (read.coordinates as GridPositions).data.buffer,
      );
      for (let i = 0; i < page.list.count; ++i) {
        expect(Array.from(GridPositions.get(page.list.coordinates, i, p))).toEqual(
          Array.from(GridPositions.get(read.coordinates, page.start + i, q)),
        );
        expect(ParticleValues.get(page.list.targets, i)).toBe(ParticleValues.get(read.targets, page.start + i));
      }
    }
    // a list without a directory: none
    expect(
      ParticleTileView.ofWhole(createParticleListFromPack(readParticlePack(writeParticlePack(list)))),
    ).toBeUndefined();
  });

  it("keeps a target's reach and volume as a reader computes them from its sites: the catalog without a pass", () => {
    const bytes = writeParticlePack(createList());
    const doc = readParticlePack(bytes).document;
    const read = createParticleListFromPack(readParticlePack(bytes));
    // the same list without them: computed from the sites
    const plain: ParticleList = {
      ...read,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    const known = ParticleTileCatalog.fromList(read).targets;
    const computed = ParticleTileCatalog.fromList(plain).targets;
    for (const id of [5, 9]) {
      const row = doc.targets.find((t) => t.id === id)!;
      expect(row.reach).toBeGreaterThan(0);
      expect(row.volume).toBeGreaterThan(0);
      expect(known.get(id)).toEqual(computed.get(id));
    }
    // a target without sites none
    expect(doc.targets.find((t) => t.id === 2)!.reach).toBeUndefined();
  });

  it('writes a tile-ordered pack: by tile in Morton order, a tile by target, the targets runs, a directory', () => {
    const list = createSpread();
    const read = createParticleListFromPack(readParticlePack(writeParticlePack(list, { tiles: 4 })));
    expect(read.count).toBe(list.count);
    const sourceOf = (p: number) => read.attributes!.get('source')!.column.value(p);
    // every particle once, its values its source's
    const seen = new Set<number>();
    for (let p = 0; p < read.count; ++p) {
      const i = sourceOf(p);
      seen.add(i);
      expect(ParticleValues.get(read.targets, p)).toBe(list.targets[i]);
      expect(ParticleValues.get(read.entities!, p)).toBe(list.entities![i]);
      expect(ParticleValues.get(read.radii!, p)).toBe(list.radii![i]);
      expectPosition(read, p, list, i);
    }
    expect(seen.size).toBe(list.count);
    // the fiber's particles at their new positions, in their order
    expect(Array.from(read.fibers!.indices, sourceOf)).toEqual(Array.from(list.fibers!.indices));
    // the targets' runs a target in a tile: none of the same target next to each other within a tile
    const directory = ParticleTiles.getDirectory(read)!;
    expect(directory.dimensions).toBe(4);
    const grid = ParticleTileGrid.create(directory.min, directory.max, [4, 4, 4]);
    const runs = read.targets as { offsets: Int32Array; values: Int32Array };
    const p3 = Vec3();
    let end = 0;
    let lastCode = -1;
    for (let k = 0; k < directory.ranges.length / 3; ++k) {
      const [id, first, count] = directory.ranges.subarray(k * 3, k * 3 + 3);
      expect(first).toBe(end);
      end = first + count;
      // in Morton order of their cells
      const code = mortonOrder3d(id % 4, Math.floor(id / 4) % 4, Math.floor(id / 16));
      expect(code).toBeGreaterThan(lastCode);
      lastCode = code;
      const box = directory.boxes.subarray(k * 6, k * 6 + 6);
      const targets: number[] = [];
      for (let p = first; p < end; ++p) {
        GridPositions.get(read.coordinates, p, p3);
        expect(ParticleTileGrid.getIdOf(grid, p3[0], p3[1], p3[2])).toBe(id);
        for (let a = 0; a < 3; ++a) {
          expect(p3[a]).toBeGreaterThanOrEqual(box[a]);
          expect(p3[a]).toBeLessThanOrEqual(box[a + 3]);
        }
        const t = ParticleValues.get(read.targets, p);
        if (targets.at(-1) !== t) targets.push(t);
      }
      // by target: each once
      expect(new Set(targets).size).toBe(targets.length);
      // its runs start at its first particle
      expect(runs.offsets).toContain(first);
    }
    expect(end).toBe(read.count);
    // no bounds stored a run (computed where asked): the targets' bounds rows none
    expect(readParticlePack(writeParticlePack(list, { tiles: 4 })).document.targets.every((t) => !t.bounds)).toBe(true);
    expect(() => writeParticlePack(list, { tiles: 3 })).toThrow(/power of two/);
  });

  it("makes a tile-ordered pack's tile sets from its directory: its tiles' grid and coarser, a tile a range", () => {
    const list = createSpread();
    const read = createParticleListFromPack(readParticlePack(writeParticlePack(list, { tiles: 4 })));
    // the same list without its directory: tiled by a pass, ordered again
    const plain: ParticleList = {
      ...read,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    expect(ParticleTiles.getDirectory(plain)).toBeUndefined();
    for (const d of [4, 2, 1, 3]) {
      const fromDirectory = ParticleTiles.fromList(read, [d, d, d]);
      const byPass = ParticleTiles.fromList(plain, [d, d, d]);
      expect([...fromDirectory.tiles.keys()]).toEqual([...byPass.tiles.keys()]);
      for (const [id, tile] of fromDirectory.tiles) {
        const other = byPass.tiles.get(id)!;
        expect(tile.count).toBe(other.count);
        expect(Array.from(tile.targets)).toEqual(Array.from(other.targets));
        expect(tile.box).toEqual(other.box);
        const a = fromDirectory.provider.getTile(id, 0) as ParticleList;
        const b = byPass.provider.getTile(id, 0) as ParticleList;
        const pa = Vec3();
        const pb = Vec3();
        for (let i = 0; i < a.count; ++i) {
          expect(Array.from(GridPositions.get(a.coordinates, i, pa))).toEqual(
            Array.from(GridPositions.get(b.coordinates, i, pb)),
          );
          expect(ParticleValues.get(a.targets, i)).toBe(ParticleValues.get(b.targets, i));
        }
        // its bounds the directory's boxes' (float32, rounded outward): the pass's, a float32 step larger at most
        const ba = Particle.getBoundary(a).box;
        const bb = Particle.getBoundary(b).box;
        for (let k = 0; k < 3; ++k) {
          expect(ba.min[k]).toBeLessThanOrEqual(bb.min[k]);
          expect(ba.max[k]).toBeGreaterThanOrEqual(bb.max[k]);
          expect(bb.min[k] - ba.min[k]).toBeLessThan(1e-4);
          expect(ba.max[k] - bb.max[k]).toBeLessThan(1e-4);
        }
        if (d === 3) continue;
        // a range of the pack: its keys its particles' indices in it, its positions a view of its block
        expect(Array.from(a.keys!)).toEqual(Array.from({ length: tile.count }, (_, i) => tile.base + i));
        const words = (a.coordinates as GridPositions).data;
        expect(words.buffer).toBe((read.coordinates as GridPositions).data.buffer);
      }
      // its sphere the same tiles' (an approximate union: the boxes a float32 step apart)
      expect(Math.abs(fromDirectory.boundingSphere!.radius / byPass.boundingSphere!.radius - 1)).toBeLessThan(1e-3);
    }
  });

  it("packs an mmCIF assembly's list: its particles and its targets' sites as read", async () => {
    const parsed = await CIF.parseText(assemblyCif).run();
    if (parsed.isError) throw new Error(parsed.message);
    const list = await createParticleListFromMmcifAssembly(RuntimeContext.Synchronous, parsed.result, {
      assemblyId: '1',
      variant: 'standard',
    });
    const read = createParticleListFromPack(readParticlePack(writeParticlePack(list)));
    const order = getPackOrder(list);
    expect(read.count).toBe(list.count);
    for (let p = 0; p < read.count; ++p) {
      const i = order[p];
      expect(ParticleValues.get(read.targets!, p)).toBe(list.targets[i]);
      expectPosition(read, p, list, i);
      // snorm16 (`ParticleRotations`)
      expect(Array.from(read.rotations!.subarray(p * 4, p * 4 + 4))).toEqual(
        Array.from(ParticleRotations.toSnorm16(list.rotations!).subarray(i * 4, i * 4 + 4)),
      );
      expect(ParticleValues.get(read.radii!, p)).toBe(list.radii![i]);
    }
    for (const [id, target] of list.targetMapping!) {
      const packed = read.targetMapping!.get(id)!;
      if (target.kind !== 'sites' || packed.kind !== 'sites') throw new Error('sites expected');
      expectClose(packed.sites.positions, target.sites.positions, packed.sites.quantized!.scale / 2 + 1e-6);
      expect(Array.from(packed.sites.radii)).toEqual(Array.from(target.sites.radii));
      // the elements as read from `type_symbol`
      expect(target.sites.elements).toBeDefined();
      expect(Array.from(packed.sites.elements!)).toEqual(Array.from(target.sites.elements!));
      expect(Array.from(packed.sites.chains!.offsets)).toEqual(Array.from(target.sites.chains!.offsets));
      expect(packed.sites.chains!.labels).toEqual(target.sites.chains!.labels);
      expect(Array.from(packed.center!)).toEqual(Array.from(ParticleTarget.getCenter(target)!));
    }
  });
});

// two chains (targets), their gen rows interleaved: the pack sorts their particles by chain
const assemblyCif = `data_pack_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA ALA A 1 1 0 0 0 1
ATOM 2 N N ALA A 1 1 1 0 0 1
ATOM 3 O O GLY B 2 1 5 5 5 1

loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 (1,2) A
1 (3) B
1 (4) A

loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 'identity operation' 1 0 0 0 0 1 0 0 0 0 1 0
2 'general operation' 0 -1 0 10 1 0 0 0 0 0 1 0
3 'general operation' 1 0 0 0 0 0 -1 5 0 1 0 0
4 'general operation' -1 0 0 -4 0 -1 0 2 0 0 1 1
`;

describe('particle pack templates', () => {
  /** four copies of `template` on a 2x2 grid (one turned 90° about z), about the origin, a leaf particle of target 3 */
  function createNested(template: ParticleList, other?: ParticleList): ParticleList {
    const turn = [0, 0, Math.SQRT1_2, Math.SQRT1_2];
    return {
      label: 'nested',
      count: 5,
      coordinates: Float32Array.of(0, 0, 0, 40, 0, 0, 0, 40, 0, 40, 40, 0, 20, 20, 20),
      rotations: Float32Array.of(0, 0, 0, 1, ...turn, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1),
      targets: Int32Array.of(0, 0, 1, 0, 3),
      targetInfo: new Map([
        [0, { name: 'copy' }],
        [1, { name: 'other copy' }],
        [3, { name: 'leaf' }],
      ]),
      targetMapping: new Map<number, ParticleTarget>([
        [0, { kind: 'particles', particles: template, center: Vec3.create(0, 0, 0) }],
        [1, { kind: 'particles', particles: other ?? template, center: Vec3.create(0, 0, 0) }],
      ]),
      radii: Float32Array.of(30, 30, 30, 30, 2),
      getParticleLabel: (i) => `copy ${i}`,
      sourceData: { kind: 'test', name: 'test', data: undefined },
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
  }

  /** each instance of `list` flattened: its target's and entity's names, radius, position and rotation, sorted */
  function getInstances(list: ParticleList) {
    const flat = ParticleNesting.flatten(list);
    const p = Vec3();
    const q = Quat();
    const out = [];
    for (let i = 0; i < flat.count; ++i) {
      GridPositions.get(flat.coordinates, i, p);
      ParticleRotations.get(flat.rotations, i, q);
      out.push({
        target: flat.targetInfo.get(ParticleValues.get(flat.targets, i))?.name,
        entity: flat.entities ? flat.entityInfo?.get(ParticleValues.get(flat.entities, i))?.name : undefined,
        radius: flat.radii ? ParticleValues.get(flat.radii, i) : undefined,
        position: Vec3.clone(p),
        rotation: Quat.clone(q),
      });
    }
    const key = (a: (typeof out)[number]) => `${a.target} ${a.position.map((v) => v.toFixed(1)).join(' ')}`;
    return out.sort((a, b) => (key(a) < key(b) ? -1 : key(a) > key(b) ? 1 : 0));
  }

  it('writes a list with templates (version 7): read back, its flattened instances the same as the written one’s', () => {
    const nested = createNested(createList());
    const bytes = writeParticlePack(nested);
    expect(new DataView(bytes.buffer).getUint32(4, true)).toBe(7);
    const read = createParticleListFromPack(readParticlePack(bytes));
    expect(read.count).toBe(5);
    expect(ParticleNesting.getExpandedCount(read)).toBe(4 * 7 + 1);
    // one template read for the two targets placing it, its particles views of the same bytes
    const a = ParticleNesting.getTemplate(read, 0)!;
    expect(ParticleNesting.getTemplate(read, 1)!.particles).toBe(a.particles);
    expect((a.particles.coordinates as GridPositions).data.buffer).toBe(
      (read.coordinates as GridPositions).data.buffer,
    );
    const expected = getInstances(nested);
    const actual = getInstances(read);
    expect(actual.map((x) => [x.target, x.entity, x.radius])).toEqual(
      expected.map((x) => [x.target, x.entity, x.radius]),
    );
    for (let i = 0; i < expected.length; ++i) {
      expect(Vec3.distance(actual[i].position, expected[i].position)).toBeLessThan(1e-2);
      expect(Math.abs(Quat.dot(actual[i].rotation, expected[i].rotation))).toBeCloseTo(1, 3);
    }
    // the template's sites as a list without templates reads them
    const sites = createParticleListFromPack(readParticlePack(writeParticlePack(createList()))).targetMapping!;
    expect(a.particles.targetMapping!.get(5)!.kind).toBe('sites');
    expect([...a.particles.targetMapping!.keys()].sort()).toEqual([...sites.keys()].sort());
  });

  it("writes a template by target in its levels' sample order and its boundary: read, neither computed again", () => {
    // 3,000 particles of three targets interleaved, spread over a box (a generator of its own: the same each run)
    let seed = 7;
    const random = () => ((seed = (seed * 1103515245 + 12345) % 2 ** 31) / 2 ** 31) * 80 - 40;
    const n = 3000;
    const targets = Int32Array.from({ length: n }, (_, i) => [4, 1, 8][i % 3]);
    const template: ParticleList = {
      label: 'spread',
      count: n,
      coordinates: Float32Array.from({ length: n * 3 }, random),
      targets,
      targetInfo: new Map([
        [1, { name: 'a' }],
        [4, { name: 'b' }],
        [8, { name: 'c' }],
      ]),
      radii: Float32Array.from(targets, (t) => t / 2),
      getParticleLabel: (i) => `${i}`,
      sourceData: { kind: 'test', name: 'test', data: undefined },
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    const pack = readParticlePack(writeParticlePack(createNested(template)));
    expect(pack.document.sampleOrdered).toBeUndefined();
    expect(pack.document.templates![0].sampleOrdered).toBe(true);
    const read = createParticleListFromPack(pack);
    const tpl = ParticleNesting.getTemplate(read, 0)!.particles;
    expect(Particle.isSampleOrdered(tpl)).toBe(true);
    expect(Particle.hasBoundary(tpl)).toBe(true);
    /** the same particles, nothing set: computed */
    const bare = (l: ParticleList): ParticleList => ({
      ...l,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    });
    // the boundary as computed, rounded outward
    const stored = Particle.getBoundary(tpl);
    const computed = Particle.getBoundary(bare(tpl));
    for (let k = 0; k < 3; ++k) {
      expect(stored.box.min[k]).toBeLessThanOrEqual(computed.box.min[k]);
      expect(stored.box.min[k]).toBeGreaterThan(computed.box.min[k] - 0.002);
      expect(stored.box.max[k]).toBeGreaterThanOrEqual(computed.box.max[k]);
      expect(stored.box.max[k]).toBeLessThan(computed.box.max[k] + 0.002);
    }
    expect(Vec3.distance(stored.sphere.center, computed.sphere.center)).toBeLessThan(1e-4);
    expect(stored.sphere.radius).toBeGreaterThanOrEqual(computed.sphere.radius);
    expect(stored.sphere.radius).toBeLessThan(computed.sphere.radius + 0.003);
    // the leaf list the levels are of takes both
    const leaves = ParticleNesting.getLeafList(read, tpl);
    expect(Particle.isSampleOrdered(leaves)).toBe(true);
    expect(Particle.getBoundary(leaves)).toBe(stored);
    // the levels' picks the sorted order's: each target's the same particles
    const byTarget = (s: Particle.SampleLevel | undefined, list: ParticleList) => {
      const out = new Map<number, number[]>();
      for (const i of s!.indices) {
        const t = ParticleValues.get(list.targets, i);
        if (!out.has(t)) out.set(t, []);
        out.get(t)!.push(i);
      }
      for (const v of out.values()) v.sort((a, b) => a - b);
      return out;
    };
    for (const level of [1, 2]) {
      const fromStored = byTarget(Particle.sampleByTarget(leaves, level), leaves);
      const fromSorted = byTarget(Particle.sampleByTarget(bare(leaves), level), leaves);
      expect(fromStored).toEqual(fromSorted);
    }
  });

  it('writes a template tile-ordered (a tile by target, its sample order kept): its pick by runs the order’s', () => {
    let seed = 11;
    const random = () => ((seed = (seed * 1103515245 + 12345) % 2 ** 31) / 2 ** 31) * 80 - 40;
    const n = 3000;
    const targets = Int32Array.from({ length: n }, (_, i) => [4, 1, 8][i % 3]);
    const template: ParticleList = {
      label: 'spread',
      count: n,
      coordinates: Float32Array.from({ length: n * 3 }, random),
      targets,
      targetInfo: new Map([
        [1, { name: 'a' }],
        [4, { name: 'b' }],
        [8, { name: 'c' }],
      ]),
      radii: Float32Array.from(targets, (t) => t / 2),
      getParticleLabel: (i) => `${i}`,
      sourceData: { kind: 'test', name: 'test', data: undefined },
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    const nested = createNested(template);
    const pack = readParticlePack(writeParticlePack(nested, { templateTiles: 4 }));
    const doc = pack.document.templates![0];
    expect(doc.tiles?.dimensions).toBe(4);
    expect(doc.sampleOrdered).toBe(true);
    const read = createParticleListFromPack(pack);
    const tpl = ParticleNesting.getTemplate(read, 0)!.particles;
    // its instances the written list's (in another order)
    const expected = getInstances(nested);
    const actual = getInstances(read);
    expect(actual.map((x) => [x.target, x.radius])).toEqual(expected.map((x) => [x.target, x.radius]));
    for (let i = 0; i < expected.length; ++i) {
      expect(Vec3.distance(actual[i].position, expected[i].position)).toBeLessThan(1e-2);
    }
    // a target a run a tile it is in: its runs apart
    expect(ParticleValues.isRuns(tpl.targets)).toBe(true);
    const runs = tpl.targets as { values: Int32Array };
    expect(runs.values.length).toBeGreaterThan(3 * 8);
    // the boundary of every particle as written
    const bare = (l: ParticleList): ParticleList => ({
      ...l,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    });
    const stored = Particle.getBoundary(tpl).sphere;
    const computed = Particle.getBoundary(bare(tpl)).sphere;
    expect(stored.radius).toBeGreaterThanOrEqual(computed.radius);
    expect(stored.radius).toBeLessThan(computed.radius + 0.003);
    // the levels' picks by its runs (the stored order: no sort) the same as a pass over it in that order
    const leaves = ParticleNesting.getLeafList(read, tpl);
    const identity = Int32Array.from({ length: leaves.count }, (_, i) => i);
    for (const level of [1, 2]) {
      const byRuns = Particle.sampleByTarget(leaves, level)!;
      const byPass = Particle.sampleByTarget(leaves, level, identity)!;
      expect(Array.from(byRuns.indices)).toEqual(Array.from(byPass.indices));
      expect(byRuns.list.scales).toEqual(byPass.list.scales);
    }
    // a target's picks spread over the template (every octant of its box), not a corner of it
    const level1 = Particle.sampleByTarget(leaves, 1)!;
    const octants = new Set<number>();
    const p = Vec3();
    for (const i of level1.indices) {
      GridPositions.get(leaves.coordinates, i, p);
      octants.add((p[0] >= 0 ? 1 : 0) + (p[1] >= 0 ? 2 : 0) + (p[2] >= 0 ? 4 : 0));
    }
    expect(octants.size).toBe(8);
  });

  it("writes a template's coarse levels: read, not sampled, the samples its load would take", () => {
    let seed = 7;
    const random = () => ((seed = (seed * 1103515245 + 12345) % 2 ** 31) / 2 ** 31) * 80 - 40;
    const n = 3000;
    const targets = Int32Array.from({ length: n }, (_, i) => [4, 1, 8][i % 3]);
    const template: ParticleList = {
      label: 'spread',
      count: n,
      coordinates: Float32Array.from({ length: n * 3 }, random),
      rotations: Float32Array.from({ length: n * 4 }, (_, k) => (k % 4 === 3 ? 1 : 0)),
      targets,
      targetInfo: new Map([
        [1, { name: 'a' }],
        [4, { name: 'b' }],
        [8, { name: 'c' }],
      ]),
      getParticleLabel: (i) => `${i}`,
      sourceData: { kind: 'test', name: 'test', data: undefined },
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    const bytes = writeParticlePack(createNested(template), { templateTiles: 4 });
    expect(readParticlePack(bytes).document.templates![0].templateLevels).toBe(ParticleNesting.MaxTemplateLevels);
    const read = () => {
      const list = createParticleListFromPack(readParticlePack(bytes));
      return { list, tpl: ParticleNesting.getTemplate(list, 0)!.particles };
    };
    const a = read();
    expect(Particle.getStoredLevels(a.tpl)?.length).toBe(ParticleNesting.MaxTemplateLevels);
    const stored = ParticleNesting.getTemplateLevels(a.list, a.tpl);
    // the same pack's levels sampled (its stored ones left out)
    const b = read();
    delete b.tpl._propertyData[Particle.StoredLevelsDescriptor.name];
    const sampled = ParticleNesting.getTemplateLevels(b.list, b.tpl);
    expect(sampled.length).toBeGreaterThan(0);
    expect(stored.map((l) => l.level)).toEqual(sampled.map((l) => l.level));
    const p = Vec3();
    const q = Vec3();
    for (let k = 0; k < stored.length; ++k) {
      const [x, y] = [stored[k].list, sampled[k].list];
      expect(Array.from(stored[k].indices)).toEqual(Array.from(sampled[k].indices));
      expect(x.targets).toEqual(y.targets);
      expect(x.scales).toEqual(y.scales);
      expect(Array.from(x.rotations!)).toEqual(Array.from(y.rotations!));
      for (let i = 0; i < x.count; ++i) {
        GridPositions.get(x.coordinates, i, p);
        GridPositions.get(y.coordinates, i, q);
        expect(Vec3.distance(p, q)).toBe(0);
      }
    }
  });

  it('warns of a template an older writer wrote (by target, its levels not stored), not of one written now', () => {
    // more than a tile's particles (`templateTilesOf`: tile-ordered by default)
    let seed = 3;
    const random = () => ((seed = (seed * 1103515245 + 12345) % 2 ** 31) / 2 ** 31) * 200 - 100;
    const n = 40_000;
    const template: ParticleList = {
      label: 'spread',
      count: n,
      coordinates: Float32Array.from({ length: n * 3 }, random),
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
    const pack = readParticlePack(writeParticlePack(createNested(template)));
    expect(pack.document.templates![0].tiles).toBeDefined();
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    try {
      createParticleListFromPack(pack);
      expect(warn).not.toHaveBeenCalled();
      // by target (tiles 1: as before tile order)
      createParticleListFromPack(readParticlePack(writeParticlePack(createNested(template), { templateTiles: 1 })));
      expect(warn).toHaveBeenCalledTimes(1);
      expect(String(warn.mock.calls[0][0])).toMatch(/not tile-ordered/);
      // its levels not stored (as before they were)
      const [t] = pack.document.templates!;
      const older = { ...t, templateLevels: undefined };
      createParticleListFromPack({ ...pack, document: { ...pack.document, templates: [older] } });
      expect(warn).toHaveBeenCalledTimes(2);
      expect(String(warn.mock.calls[1][0])).toMatch(/coarse levels not stored/);
    } finally {
      warn.mockRestore();
    }
  });

  it('writes each template once, a list without templates as version 6', () => {
    const template = createList();
    const once = writeParticlePack(createNested(template));
    const twice = writeParticlePack(createNested(template, createList()));
    const templates = (b: Uint8Array) => readParticlePack(b).document.templates?.length;
    expect(templates(once)).toBe(1);
    expect(templates(twice)).toBe(2);
    expect(twice.byteLength).toBeGreaterThan(once.byteLength);
    const plain = writeParticlePack(template);
    expect(new DataView(plain.buffer).getUint32(4, true)).toBe(6);
    expect(readParticlePack(plain).document.templates).toBeUndefined();
  });

  it('refuses a list with templates tile-ordered, nesting deeper than a level, a template not in the pack', () => {
    const nested = createNested(createList());
    expect(() => writeParticlePack(nested, { tiles: 4 })).toThrow(/not tile-ordered/);
    expect(() => writeParticlePack(createNested(nested))).toThrow(/templates/);
    const doc = (targets: unknown[], templates?: unknown[]) =>
      new TextEncoder().encode(
        JSON.stringify({
          count: 0,
          entities: [],
          compartments: [],
          attributes: [],
          siteRadii: [],
          blocks: [],
          targets,
          templates,
        }),
      );
    const target = { id: 0, entity: -1, compartment: -1, particleOffset: 0, particleCount: 0 };
    const placing = { ...target, template: { index: 1, center: [0, 0, 0] } };
    expect(() => parseParticlePackDocument(doc([placing]), 0)).toThrow(/has none/);
    const template = JSON.parse(new TextDecoder().decode(doc([])));
    expect(() => parseParticlePackDocument(doc([placing], [template]), 0)).toThrow(/not a template/);
    const deeper = { ...template, templates: [template] };
    expect(() => parseParticlePackDocument(doc([], [deeper]), 0)).toThrow(/templates of its own/);
  });
});
