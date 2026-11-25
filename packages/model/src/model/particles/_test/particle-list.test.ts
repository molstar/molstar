/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { CustomProperties } from '../../custom-property.js';
import { OrderedSet } from '@molstar/core/data/int';
import { lociLabel } from '../../label.js';
import {
  getParticleSpacingRadius,
  getParticleTargetGroups,
  getParticleTransforms,
  ParticleRotations,
  Particle,
  ParticleList,
  type ParticleRuns,
  ParticleSites,
  ParticleTarget,
  ParticleValues,
} from '../particle-list.js';
import { Grid } from '../../volume/grid.js';
import type { Volume } from '../../volume/volume.js';
import { grid } from '../../volume/_test/grid-test-util.js';

function createParticles(coordinates: Float32Array, rotations?: ParticleRotations): ParticleList {
  const count = coordinates.length / 3;
  return {
    count,
    keys: Int32Array.from({ length: count }, (_, i) => i),
    targets: new Int32Array(count),
    targetInfo: new Map([[0, {}]]),
    coordinates,
    rotations,
    getParticleLabel: (index) => `${index}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

describe('ParticleTarget.getCenter', () => {
  it("places a volume by its box center, whatever its isosurface's", () => {
    // density only in one corner: the isosurface's sphere is far from the box's center
    const g = grid([0, 1, 2], Float32Array, [4, 6, 8], (i, j, k) => (i < 1 && j < 1 && k < 1 ? 1 : 0));
    const matrix = Mat4.mul(
      Mat4(),
      Mat4.fromTranslation(Mat4(), Vec3.create(1, 2, 3)),
      Mat4.fromScaling(Mat4(), Vec3.create(2, 2, 2)),
    );
    const volume = { grid: { ...g, transform: { kind: 'matrix', matrix } } } as Volume;
    const center = ParticleTarget.getCenter({ kind: 'volume', volume })!;
    expect(Array.from(center)).toEqual([1 + 4, 2 + 6, 3 + 8]);
    const iso = Grid.getIsosurfaceBoundingSphere(volume.grid, 0.5);
    expect(Vec3.distance(iso.center, center)).toBeGreaterThan(5);
  });
});

describe('getParticleTransforms', () => {
  it('uses identity bases without rotations', () => {
    const particles = createParticles(new Float32Array([1, 2, 3, 4, 5, 6]));
    const t = getParticleTransforms(particles);

    const expected = Mat4.fromTranslation(Mat4(), Vec3.create(1, 2, 3));
    for (let j = 0; j < 16; ++j) expect(t[j]).toBeCloseTo(expected[j], 6);
    expect(t[16 + 12]).toBe(4);
    expect(t[16 + 13]).toBe(5);
    expect(t[16 + 14]).toBe(6);
  });

  it('matches Mat4.fromQuat for rotated particles', () => {
    // deterministic LCG so the test is reproducible
    let seed = 3;
    const rand = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;

    const count = 25;
    const coordinates = new Float32Array(count * 3);
    const rotations = new Float32Array(count * 4);
    const q = Quat();
    for (let i = 0; i < count; ++i) {
      coordinates[i * 3] = rand() * 100;
      coordinates[i * 3 + 1] = rand() * 100;
      coordinates[i * 3 + 2] = rand() * 100;
      Quat.set(q, rand() - 0.5, rand() - 0.5, rand() - 0.5, rand() - 0.5);
      Quat.normalize(q, q);
      rotations.set(q, i * 4);
    }

    const particles = createParticles(coordinates, rotations);
    const t = getParticleTransforms(particles);

    const m = Mat4();
    for (let i = 0; i < count; ++i) {
      Quat.set(q, rotations[i * 4], rotations[i * 4 + 1], rotations[i * 4 + 2], rotations[i * 4 + 3]);
      Mat4.fromQuat(m, q);
      m[12] = coordinates[i * 3];
      m[13] = coordinates[i * 3 + 1];
      m[14] = coordinates[i * 3 + 2];
      for (let j = 0; j < 16; ++j) {
        expect(t[i * 16 + j]).toBeCloseTo(m[j], 5);
      }
    }
  });

  it('caches the result on the particle list', () => {
    const particles = createParticles(new Float32Array([1, 2, 3]));
    expect(getParticleTransforms(particles)).toBe(getParticleTransforms(particles));
  });
});

describe('getParticleSpacingRadius', () => {
  it('is half the distance to the nearest particle, for particles on a lattice', () => {
    const coords: number[] = [];
    for (let i = 0; i < 6; i++)
      for (let j = 0; j < 6; j++) for (let k = 0; k < 6; k++) coords.push(i * 10, j * 10, k * 10);
    const list = createParticles(new Float32Array(coords));
    expect(getParticleSpacingRadius(list)).toBeCloseTo(5, 4);
    // cached
    expect(getParticleSpacingRadius(list)).toBe(getParticleSpacingRadius(list));
  });

  it('is 1 for a single particle or particles all at one point', () => {
    expect(getParticleSpacingRadius(createParticles(new Float32Array([1, 2, 3])))).toBe(1);
    expect(getParticleSpacingRadius(createParticles(new Float32Array([1, 2, 3, 1, 2, 3])))).toBe(1);
  });
});

describe('ParticleRotations', () => {
  // deterministic LCG so the test is reproducible
  let seed = 7;
  const rand = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
  const count = 200;
  const rotations = new Float32Array(count * 4);
  const q = Quat();
  for (let i = 0; i < count; ++i) {
    Quat.normalize(q, Quat.set(q, rand() - 0.5, rand() - 0.5, rand() - 0.5, rand() - 0.5));
    rotations.set(q, i * 4);
  }
  rotations.set([0, 0, 0, -1], 0);

  it('quantizes unit quaternions within 1e-4 rad, w not below 0, the identity exact', () => {
    const snorm = ParticleRotations.toSnorm16(rotations);
    expect(snorm).toBeInstanceOf(Int16Array);
    expect(ParticleRotations.toSnorm16(snorm)).toBe(snorm);
    expect(Array.from(snorm.subarray(0, 4))).toEqual([0, 0, 0, ParticleRotations.Snorm16One]);
    const a = Quat();
    const b = Quat();
    for (let i = 0; i < count; ++i) {
      expect(snorm[i * 4 + 3]).toBeGreaterThanOrEqual(0);
      ParticleRotations.get(snorm, i, a);
      ParticleRotations.get(rotations, i, b);
      expect(Quat.magnitude(a)).toBeCloseTo(1, 6);
      // the angle between them: 2 acos(|a · b|)
      expect(2 * Math.acos(Math.min(1, Math.abs(Quat.dot(a, b))))).toBeLessThan(1e-4);
    }
  });

  it('gives the transforms of the quaternions they quantize, read as they are', () => {
    const coordinates = new Float32Array(count * 3).map((_, k) => k * 0.5);
    const snorm = ParticleRotations.toSnorm16(rotations);
    const unit = new Float32Array(count * 4);
    for (let i = 0; i < count; ++i) unit.set(ParticleRotations.get(snorm, i, Quat()), i * 4);
    const a = getParticleTransforms(createParticles(coordinates, snorm));
    const b = getParticleTransforms(createParticles(coordinates, unit));
    for (let k = 0; k < a.length; ++k) expect(a[k]).toBeCloseTo(b[k], 5);
    // a move: float32 unit quaternions, those of the turned quantized ones
    const turn = Mat4.fromRotation(Mat4(), 0.4, Vec3.create(1, 1, 0));
    const ta = Particle.transform(createParticles(coordinates, snorm), turn).rotations!;
    const tb = Particle.transform(createParticles(coordinates, unit), turn).rotations!;
    expect(ta).toBeInstanceOf(Float32Array);
    for (let k = 0; k < ta.length; ++k) expect(ta[k]).toBeCloseTo(tb[k], 5);
  });
});

describe('Particle.withPoses', () => {
  it('is the same particles posed otherwise: one line of them, what is cached by position its own', () => {
    const list = createParticles(new Float32Array([1, 0, 0, 0, 2, 0, 5, 5, 5]), new Float32Array(12).fill(0.5));
    const groups = getParticleTargetGroups(list);
    const spacing = getParticleSpacingRadius(list);
    const boundary = Particle.getBoundary(list);
    const frame = Particle.withPoses(list, new Float32Array([2, 0, 0, 0, 3, 0, 5, 6, 5]));
    const next = Particle.withPoses(frame, new Float32Array([3, 0, 0, 0, 4, 0, 5, 7, 5]), new Float32Array(12));
    // one line: its first the list posed first, either of it the same particles
    expect(frame.posesFrom).toBe(list);
    expect(next.posesFrom).toBe(list);
    expect(Particle.isSamePoses(list, next)).toBe(true);
    expect(Particle.isSamePoses(frame, next)).toBe(true);
    expect(Particle.isSamePoses(list, createParticles(new Float32Array(9)))).toBe(false);
    // the poses its own, the rest the list's
    expect(frame.coordinates[0]).toBe(2);
    expect(frame.rotations).toBe(list.rotations);
    expect(next.rotations![0]).toBe(0);
    expect(frame.targets).toBe(list.targets);
    // not by position kept, by position (the spacing too: the distances change) computed again; not a move
    expect(getParticleTargetGroups(frame)).toBe(groups);
    expect(getParticleSpacingRadius(frame)).not.toBe(spacing);
    expect(Particle.getBoundary(frame)).not.toBe(boundary);
    expect(Particle.getMove(list, frame)).toBeUndefined();
    expect(() => Particle.withPoses(list, new Float32Array(6))).toThrow();
  });
});

describe('Particle.transform', () => {
  const turn = Quat.setAxisAngle(Quat(), Vec3.create(0, 0, 1), Math.PI / 2);
  const matrix = Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), Vec3.create(10, 0, 0)), Mat4.fromQuat(Mat4(), turn));

  it('moves the coordinates and gives orientations where it turns a list without them', () => {
    const list = createParticles(new Float32Array([1, 0, 0, 0, 2, 0]));
    Particle.getBoundary(list);
    const moved = Particle.transform(list, matrix);
    const c = Array.from(moved.coordinates);
    [10, 1, 0, 8, 0, 0].forEach((x, i) => expect(c[i]).toBeCloseTo(x, 5));
    expect(Array.from(moved.rotations!.subarray(0, 4))).toEqual(Array.from(turn).map(Math.fround));
    // the source untouched, the caches not shared
    expect(Array.from(list.coordinates)).toEqual([1, 0, 0, 0, 2, 0]);
    expect(Particle.getBoundary(moved)).not.toBe(Particle.getBoundary(list));
    expect(Particle.getBoundary(moved).sphere.center[0]).toBeGreaterThan(8);
  });

  it('turns the orientations a list has; a move alone keeps a list without them so', () => {
    const own = Quat.setAxisAngle(Quat(), Vec3.create(1, 0, 0), Math.PI / 3);
    const list = createParticles(new Float32Array([0, 0, 0]), Float32Array.from(own));
    const q = Particle.transform(list, matrix).rotations!;
    const expected = Quat.multiply(Quat(), turn, own);
    for (let j = 0; j < 4; j++) expect(q[j]).toBeCloseTo(expected[j], 6);

    const plain = createParticles(new Float32Array([0, 0, 0]));
    expect(Particle.transform(plain, Mat4.fromTranslation(Mat4(), Vec3.create(1, 2, 3))).rotations).toBeUndefined();
  });

  it('records the move, so two lists of a line of moves give the move between them', () => {
    const list = createParticles(new Float32Array([1, 0, 0, 0, 2, 0]));
    const second = Mat4.fromTranslation(Mat4(), Vec3.create(0, 5, 0));
    const once = Particle.transform(list, matrix);
    // moved again from the first moved: one line, the moves composed
    const twice = Particle.transform(once, second);
    expect(once.movedFrom!.particles).toBe(list);
    expect(twice.movedFrom!.particles).toBe(list);
    expect(Mat4.areEqual(Particle.getMove(list, once)!, matrix, 1e-6)).toBe(true);
    expect(Mat4.areEqual(Particle.getMove(once, twice)!, second, 1e-5)).toBe(true);
    // moving `once` by the move gives `twice`
    const p = Vec3.transformMat4(Vec3(), Vec3.fromArray(Vec3(), once.coordinates, 3), Particle.getMove(once, twice)!);
    for (let k = 0; k < 3; k++) expect(p[k]).toBeCloseTo(twice.coordinates[3 + k], 4);
    // targets attached keep the line: a list moved under them gives its move
    const withTargets = Particle.withTargets(list, new Map());
    expect(Mat4.areEqual(Particle.getMove(withTargets, Particle.withTargets(once, new Map()))!, matrix, 1e-6)).toBe(
      true,
    );
    // another list: none
    expect(Particle.getMove(createParticles(new Float32Array([1, 0, 0, 0, 2, 0])), once)).toBeUndefined();
  });

  it('keeps what it cached not by position: the target groups, the spacing of a move without scale', () => {
    const list = createParticles(new Float32Array([1, 0, 0, 0, 2, 0, 5, 5, 5]));
    const groups = getParticleTargetGroups(list);
    const spacing = getParticleSpacingRadius(list);
    const transforms = getParticleTransforms(list);
    const moved = Particle.transform(list, matrix);
    expect(getParticleTargetGroups(moved)).toBe(groups);
    expect(getParticleSpacingRadius(moved)).toBe(spacing);
    // by position: computed again
    expect(getParticleTransforms(moved)).not.toBe(transforms);
    expect(getParticleTransforms(moved)[12]).toBeCloseTo(moved.coordinates[0], 5);
    // a scale changes the distances
    const scaled = Particle.transform(list, Mat4.fromUniformScaling(Mat4(), 2));
    expect(getParticleTargetGroups(scaled)).toBe(groups);
    expect(getParticleSpacingRadius(scaled)).toBeCloseTo(spacing * 2, 5);
  });
});

describe('a particle loci label', () => {
  it('says the particle in bold, its list beneath, the list name escaped', () => {
    const list = { ...createParticles(new Float32Array([0, 0, 0, 5, 0, 0])), label: 'picks <1>.em' };
    const one = lociLabel(Particle.Loci(list, OrderedSet.ofSingleton(1)));
    expect(one).toBe('<b>1</b></br><small>picks &lt;1&gt;.em</small>');
    expect(lociLabel(Particle.Loci(list, OrderedSet.ofBounds(0, 2)))).toBe(
      '<b>2 particles</b></br><small>picks &lt;1&gt;.em</small>',
    );
  });

  it('is plain text when asked (no tags, nothing escaped)', () => {
    const list = { ...createParticles(new Float32Array([0, 0, 0, 5, 0, 0])), label: 'picks <1>.em' };
    expect(lociLabel(Particle.Loci(list, OrderedSet.ofSingleton(1)), { htmlStyling: false })).toBe('1 | picks <1>.em');
  });
});

describe('ParticleValues', () => {
  // runs: [0, 2) 7, [2, 2) empty, [2, 5) 3, [5, 6) 7
  const runs = ParticleValues.ofRuns(Int32Array.from([0, 2, 2, 5, 6]), Int32Array.from([7, 9, 3, 7]));
  const array = Int32Array.from([7, 7, 3, 3, 3, 7]);

  it('reads a particle of an array or of runs alike', () => {
    for (let i = 0; i < 6; ++i) expect(ParticleValues.get(runs, i)).toBe(array[i]);
    // in any order (the run found last tried first)
    for (const i of [5, 0, 3, 1, 4, 2]) expect(ParticleValues.get(runs, i)).toBe(array[i]);
    expect(Array.from(ParticleValues.toArray(runs, 6, Int32Array))).toEqual(Array.from(array));
    expect(ParticleValues.toArray(array, 6, Int32Array)).toBe(array);
  });

  it('visits the runs: those of runs (none empty), those of equal values of an array', () => {
    const of = (v: ParticleValues<Int32Array>) => {
      const out: number[][] = [];
      ParticleValues.forEachRun(v, 6, (start, end, value) => out.push([start, end, value]));
      return out;
    };
    const expected = [
      [0, 2, 7],
      [2, 5, 3],
      [5, 6, 7],
    ];
    expect(of(runs)).toEqual(expected);
    expect(of(array)).toEqual(expected);
  });

  it('visits the runs of two together, a run ending where either does', () => {
    const other = ParticleValues.ofRuns(Int32Array.from([0, 3, 6]), Int32Array.from([1, 2]));
    const out: number[][] = [];
    ParticleValues.forEachRun2(runs, other, 6, (start, end, a, b) => out.push([start, end, a, b]));
    expect(out).toEqual([
      [0, 2, 7, 1],
      [2, 3, 3, 1],
      [3, 5, 3, 2],
      [5, 6, 7, 2],
    ]);
    const arrays: number[][] = [];
    ParticleValues.forEachRun2(array, Int32Array.from([1, 1, 1, 2, 2, 2]), 6, (start, end, a, b) =>
      arrays.push([start, end, a, b]),
    );
    expect(arrays).toEqual(out);
  });

  it('maps a value a run for runs, a value a particle for an array', () => {
    const mapped = ParticleValues.map(runs, 6, (v) => v * 10);
    expect(ParticleValues.isRuns(mapped)).toBe(true);
    expect(Array.from(ParticleValues.toArray(mapped, 6, Int32Array))).toEqual([70, 70, 30, 30, 30, 70]);
    expect(Array.from(ParticleValues.map(array, 6, (v) => v * 10) as Int32Array)).toEqual([70, 70, 30, 30, 30, 70]);
  });

  it('maps runs of one value side by side once, each run kept', () => {
    const same = ParticleValues.ofRuns(Int32Array.of(0, 1, 3, 4, 6), Int32Array.of(5, 5, 2, 5));
    const f = vi.fn((v: number) => v + 1);
    const mapped = ParticleValues.map(same, 6, f);
    expect(f).toHaveBeenCalledTimes(3);
    expect(ParticleValues.isRuns(mapped) && Array.from(mapped.offsets)).toEqual([0, 1, 3, 4, 6]);
    expect(Array.from(ParticleValues.toArray(mapped, 6, Int32Array))).toEqual([6, 6, 6, 3, 6, 6]);
  });

  it('groups targets of runs as of an array: a target of one run an interval', () => {
    const base = createParticles(new Float32Array(18));
    const fromRuns = getParticleTargetGroups({ ...base, targets: runs, _propertyData: Object.create(null) });
    const fromArray = getParticleTargetGroups({ ...base, targets: array, _propertyData: Object.create(null) });
    expect(Array.from(fromRuns.targetIds)).toEqual([3, 7]);
    expect(Array.from(fromRuns.targetIds)).toEqual(Array.from(fromArray.targetIds));
    for (let g = 0; g < 2; ++g) {
      expect(OrderedSet.toArray(fromRuns.sets[g])).toEqual(OrderedSet.toArray(fromArray.sets[g]));
    }
    // target 3 one run: an interval, no index a particle
    expect(OrderedSet.isInterval(fromRuns.sets[0])).toBe(true);
    expect(Particle.is({ ...base, targets: runs, keys: undefined })).toBe(true);
  });
});

describe('a target in several runs (a tile-ordered list: a run a tile it is in)', () => {
  // 12 particles in four tiles of three: targets 1 and 2 alternating by tile, runs [0,3) 1, [3,6) 2, [6,9) 1, [9,12) 2
  const coordinates = new Float32Array(36).map((_, k) => ((k * 7) % 11) - 5);
  const offsets = Int32Array.of(0, 3, 6, 9, 12);
  const list: ParticleList = {
    ...createParticles(coordinates),
    targets: ParticleValues.ofRuns(offsets, Int32Array.of(1, 2, 1, 2)),
    targetInfo: new Map([
      [1, {}],
      [2, {}],
    ]),
  };
  const one = OrderedSet.ofSortedArray(Int32Array.of(0, 1, 2, 6, 7, 8));

  it('walks the runs of its particles alone, with their counts', () => {
    const seen: [number, number][] = [];
    ParticleValues.forEachRunOf(list.targets as ParticleRuns, one, (r, c) => {
      seen.push([r, c]);
    });
    expect(seen).toEqual([
      [0, 3],
      [2, 3],
    ]);
    // a part of a run: its count of it; stopped where asked
    seen.length = 0;
    ParticleValues.forEachRunOf(
      list.targets as ParticleRuns,
      OrderedSet.ofSortedArray(Int32Array.of(1, 4, 5, 10)),
      (r, c) => {
        seen.push([r, c]);
        return r < 1;
      },
    );
    expect(seen).toEqual([
      [0, 1],
      [1, 2],
    ]);
  });

  it("bounds them by their runs' bounds united, none where they are not whole runs", () => {
    const bounds = Particle.getIndicesBounds(list, one)!;
    const p = Vec3();
    const min = Vec3.create(Infinity, Infinity, Infinity);
    const max = Vec3.create(-Infinity, -Infinity, -Infinity);
    for (let k = 0; k < OrderedSet.size(one); ++k) {
      Vec3.fromArray(p, coordinates, OrderedSet.getAt(one, k) * 3);
      Vec3.min(min, min, p);
      Vec3.max(max, max, p);
    }
    expect(Array.from(bounds.min)).toEqual(Array.from(min));
    expect(Array.from(bounds.max)).toEqual(Array.from(max));
    // every particle within the distance of the box's center
    const c = Vec3.scale(Vec3(), Vec3.add(Vec3(), min, max), 0.5);
    for (let k = 0; k < OrderedSet.size(one); ++k) {
      Vec3.fromArray(p, coordinates, OrderedSet.getAt(one, k) * 3);
      expect(Vec3.distance(p, c)).toBeLessThanOrEqual(bounds.distance + 1e-9);
    }
    // one run: its own, as getRangeBounds
    expect(Particle.getIndicesBounds(list, OrderedSet.ofRange(3, 5))).toEqual(Particle.getRangeBounds(list, 3, 6));
    expect(Particle.getIndicesBounds(list, OrderedSet.ofSortedArray(Int32Array.of(0, 1, 6, 7, 8)))).toBeUndefined();
  });
});

describe('ParticleSites.coarsen', () => {
  it("keeps every stride-th site, a chain's first its first kept one, the strides multiplied", () => {
    // 7 sites in chains of 3, 1 and 3: x the site's index
    const data = new Int16Array(7 * 4);
    for (let i = 0; i < 7; ++i) data.set([i, 2 * i, 3 * i, i % 2], i * 4);
    const chains = { offsets: Int32Array.of(0, 3, 4, 7), labels: ['A', 'B', 'C'] };
    const extremes = Int16Array.of(0, 0, 0, 2, 4, 6, 3, 6, 9, 3, 6, 9, 4, 8, 12, 6, 12, 18);
    const typeRadii = Float32Array.of(1, 2);
    const sites = ParticleSites.ofQuantized(
      7,
      { data, origin: Vec3.create(1, 0, 0), scale: 0.5, typeRadii, extremes },
      'x',
      chains,
    );

    const coarse = ParticleSites.coarsen(sites, 2);
    expect(coarse.count).toBe(4);
    expect(coarse.stride).toBe(2);
    expect(Array.from({ length: 4 }, (_, i) => ParticleSites.getPosition(coarse, i, Vec3())[0])).toEqual([1, 2, 3, 4]);
    expect(Array.from({ length: 4 }, (_, i) => ParticleSites.getRadius(coarse, i))).toEqual([1, 1, 1, 1]);
    // chains: sites 0, 2 | (3 dropped) | 4, 6
    expect(Array.from(coarse.chains!.offsets)).toEqual([0, 2, 2, 4]);
    expect(coarse.quantized!.extremes).toBe(extremes);
    expect(coarse.label).toBe('x');
    // nested: every 4th of the whole
    expect(ParticleSites.coarsen(coarse, 2).stride).toBe(4);
    expect(ParticleSites.coarsen(sites, 1)).toBe(sites);

    const plain = ParticleSites.coarsen(
      {
        count: 5,
        positions: Float32Array.from({ length: 15 }, (_, k) => k),
        radii: Float32Array.of(1, 2, 3, 4, 5),
        elements: Uint8Array.of(6, 7, 8, 1, 6),
      },
      3,
    );
    expect([plain.count, plain.stride, Array.from(plain.radii), Array.from(plain.elements!)]).toEqual([
      2,
      3,
      [1, 4],
      [6, 1],
    ]);
    expect(Array.from(plain.positions)).toEqual([0, 1, 2, 9, 10, 11]);
  });
});
