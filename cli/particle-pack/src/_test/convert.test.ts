/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { CustomProperties } from '@molstar/model/model/custom-property';
import {
  type ParticleList,
  ParticleRotations,
  ParticleSites,
  type ParticleTarget,
  ParticleValues,
} from '@molstar/model/model/particles/particle-list';
import { VdwRadius } from '@molstar/model/model/structure/model/properties/atomic/measures';
import type { ElementSymbol } from '@molstar/model/model/structure/model/types';
import { createParticleListFromPack, readParticlePack, writeParticlePack } from '@molstar/model/formats/particles/pack';
import { mmcifToPack, packToMmcif, packToPack, parseMmcif } from '../convert.js';
import { ParticleTiles } from '@molstar/model/model/particles/particle-tiles';

const vdw = (s: string) => Math.fround(VdwRadius(s as unknown as ElementSymbol));

/**
 * A list as a pack has it (targets as ranges): `cellpack`, three single-chain targets with compartments, one of them a
 * fiber's (two fibers); else two targets of two chains each (PetWorld's models, named as their entities)
 */
function createList(cellpack: boolean): ParticleList {
  const counts = cellpack ? [3, 4, 2] : [3, 2];
  const n = counts.reduce((a, b) => a + b, 0);
  const targets = new Int32Array(n);
  const coordinates = new Float32Array(n * 3);
  const rotations = new Float32Array(n * 4);
  const q = Quat();
  let p = 0;
  counts.forEach((c, t) => {
    for (let k = 0; k < c; ++k, ++p) {
      targets[p] = t;
      coordinates.set([p * 10, -p * 3, p * 0.5], p * 3);
      Quat.setAxisAngle(q, Vec3.normalize(Vec3(), Vec3.create(1, p, 2)), p * 0.7);
      Quat.toArray(q, rotations, p * 4);
    }
  });
  const targetMapping = new Map<number, ParticleTarget>();
  // sulfur and phosphorus share a radius: told apart by their elements
  const symbols = ['C', 'N', 'O', 'S', 'P'];
  const numbers = [6, 7, 8, 16, 15];
  counts.forEach((_, t) => {
    const count = 3 + t;
    const positions = new Float32Array(count * 3).map((_, k) => k * 0.75 + t);
    const radii = Float32Array.from({ length: count }, (_, k) => vdw(symbols[(k + t) % symbols.length]));
    const elements = Uint8Array.from({ length: count }, (_, k) => numbers[(k + t) % numbers.length]);
    const chains = cellpack
      ? { offsets: Int32Array.from([0, count]), labels: [`C${t}`] }
      : { offsets: Int32Array.from([0, 2, count]), labels: ['A', 'B'] };
    const sites = { count, positions, radii, elements, label: `t${t}`, chains };
    // the center the reader computes: the sites' centroid
    const center = Vec3();
    for (let k = 0; k < count; ++k) Vec3.add(center, center, Vec3.fromArray(Vec3(), positions, k * 3));
    targetMapping.set(t, { kind: 'sites', sites, center: Vec3.scale(center, center, 1 / count) });
  });
  return {
    entryId: 'TEST',
    count: n,
    keys: Int32Array.from({ length: n }, (_, i) => i),
    targets,
    targetInfo: new Map(counts.map((_, t) => [t, { name: cellpack ? `C${t}` : `entity${t}`, entity: t }])),
    targetMapping,
    entities: Int32Array.from(targets),
    entityInfo: new Map(counts.map((_, t) => [t, { name: `entity${t}` }])),
    compartments: cellpack ? Int32Array.from(targets, (t) => (t === 1 ? 1 : 0)) : undefined,
    compartmentInfo: cellpack
      ? new Map([
          [0, { name: 'root.cyto' }],
          [1, { name: 'root.cyto.fibers' }],
        ])
      : undefined,
    coordinates,
    rotations,
    // target 1's particles as two fibers of two
    fibers: cellpack
      ? { count: 2, offsets: Int32Array.from([0, 2, 4]), indices: Int32Array.from([3, 4, 5, 6]) }
      : undefined,
    getParticleLabel: (i) => `p${i}`,
    sourceData: { kind: 'test', name: 'test', data: undefined },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

/** `actual` within `eps` of `expected`, value by value (quantized sites: within half the quantization step) */
function expectClose(actual: ArrayLike<number>, expected: ArrayLike<number>, eps: number) {
  expect(actual.length).toBe(expected.length);
  for (let k = 0; k < expected.length; ++k) expect(Math.abs(actual[k] - expected[k])).toBeLessThanOrEqual(eps);
}

const pa = Vec3();
const pb = Vec3();
const qa = Quat();
const qb = Quat();

/** `b` (a pack's list read from its mmCIF) as `a`: the particles in order, their poses, the targets' sites */
function expectSameList(a: ParticleList, b: ParticleList) {
  expect(b.count).toBe(a.count);
  for (let p = 0; p < a.count; ++p) {
    GridPositions.get(a.coordinates, p, pa);
    GridPositions.get(b.coordinates, p, pb);
    for (let k = 0; k < 3; ++k) expect(pb[k]).toBeCloseTo(pa[k], 3);
    // a quaternion up to its sign
    const dot = Quat.dot(ParticleRotations.get(a.rotations, p, qa), ParticleRotations.get(b.rotations, p, qb));
    expect(Math.abs(dot)).toBeCloseTo(1, 5);
    expect(b.targetInfo.get(ParticleValues.get(b.targets!, p))!.name).toBe(
      a.targetInfo.get(ParticleValues.get(a.targets!, p))!.name,
    );
    expect(b.entityInfo!.get(ParticleValues.get(b.entities!, p))!.name).toBe(
      a.entityInfo!.get(ParticleValues.get(a.entities!, p))!.name,
    );
    if (a.compartments) {
      expect(b.compartmentInfo!.get(ParticleValues.get(b.compartments!, p))!.name).toBe(
        a.compartmentInfo!.get(ParticleValues.get(a.compartments!, p))!.name,
      );
    }
  }
  for (const [id, target] of a.targetMapping!) {
    const name = a.targetInfo.get(id)!.name;
    const bid = [...b.targetInfo].find(([, info]) => info.name === name)?.[0];
    const other = b.targetMapping!.get(bid!)!;
    if (target.kind !== 'sites' || other.kind !== 'sites') throw new Error('sites expected');
    expectClose(other.sites.positions, target.sites.positions, target.sites.quantized!.scale / 2 + 1e-6);
    expect(Array.from(other.sites.radii)).toEqual(Array.from(target.sites.radii));
    const elementsOf = (s: ParticleSites) => Array.from({ length: s.count }, (_, k) => ParticleSites.getElement(s, k));
    expect(elementsOf(other.sites)).toEqual(elementsOf(target.sites));
    expect(Array.from(other.sites.chains!.offsets)).toEqual(Array.from(target.sites.chains!.offsets));
    // the centroid of the sites as read back (quantized), the center the pack keeps
    for (let k = 0; k < 3; ++k) expect(other.center![k]).toBeCloseTo(target.center![k], 3);
  }
  if (a.fibers) {
    expect(b.fibers!.count).toBe(a.fibers.count);
    expect(Array.from(b.fibers!.offsets)).toEqual(Array.from(a.fibers.offsets));
    expect(Array.from(b.fibers!.indices)).toEqual(Array.from(a.fibers.indices));
  }
}

describe('particle pack converter', () => {
  for (const cellpack of [true, false]) {
    for (const binary of [true, false]) {
      it(`writes a pack as ${binary ? 'binary' : 'text'} mmCIF the particle reader reads as its list (${cellpack ? "cellPACK's layout" : "PetWorld's"})`, async () => {
        const pack = writeParticlePack(createList(cellpack));
        const list = createParticleListFromPack(readParticlePack(pack));
        const { data, unmatchedRadii } = packToMmcif(pack, binary);
        expect(unmatchedRadii).toBe(0);
        const bytes = typeof data === 'string' ? new TextEncoder().encode(data) : data;
        const back = createParticleListFromPack(readParticlePack(await mmcifToPack(await parseMmcif(bytes, binary))));
        expectSameList(list, back);
      });
    }
  }

  it('counts the sites of a radius no element has, their elements not known (written as carbon)', () => {
    const list = createList(true);
    const target = list.targetMapping!.get(0)!;
    if (target.kind !== 'sites') throw new Error('sites expected');
    target.sites.radii[0] = 0.123;
    list.targetMapping!.set(0, { ...target, sites: { ...target.sites, elements: undefined } });

    expect(packToMmcif(writeParticlePack(list), true).unmatchedRadii).toBe(1);
  });

  it('writes a tile-ordered pack (with levels) as mmCIF: its particles by target, the same set', async () => {
    const list = createList(true);
    const pack = writeParticlePack(list, { tiles: 2, levels: true });
    const { data } = packToMmcif(pack, true);
    const back = createParticleListFromPack(
      readParticlePack(await mmcifToPack(await parseMmcif(data as Uint8Array, true))),
    );
    expect(back.count).toBe(list.count);
    const key = (l: ParticleList) =>
      Array.from({ length: l.count }, (_, i) =>
        Array.from(GridPositions.get(l.coordinates, i, Vec3()), (v) => v.toFixed(1)).join(','),
      ).sort();
    expect(key(back)).toEqual(key(createParticleListFromPack(readParticlePack(writeParticlePack(list)))));
  });

  it('writes a pack again, tile-ordered with tiles: the same particles, a tile directory', () => {
    const pack = writeParticlePack(createList(true));
    const tiled = createParticleListFromPack(readParticlePack(packToPack(pack, { tiles: 2 })));
    const plain = createParticleListFromPack(readParticlePack(packToPack(pack)));
    expect(tiled.count).toBe(plain.count);
    expect(ParticleTiles.getDirectory(tiled)?.dimensions).toBe(2);
    expect(ParticleTiles.getDirectory(plain)).toBeUndefined();
    // the same particles in another order: the same set of positions
    const key = (l: typeof tiled) =>
      Array.from({ length: l.count }, (_, i) =>
        Array.from(GridPositions.get(l.coordinates, i, Vec3())).join(','),
      ).sort();
    expect(key(tiled)).toEqual(key(plain));
  });
});
