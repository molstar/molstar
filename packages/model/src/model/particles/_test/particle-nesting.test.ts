/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { CustomProperties } from '../../custom-property.js';
import {
  getParticleTransforms,
  Particle,
  type ParticleList,
  ParticleRotations,
  type ParticleTarget,
  ParticleValues,
} from '../particle-list.js';
import { ParticleNesting, type ParticlePath } from '../particle-nesting.js';
import { ParticleDensityGrid } from '../particle-tile-opaque.js';

/** a leaf target (never drawn here) */
const leaf = (name: string): ParticleTarget => ({ kind: 'shape', shape: { name } as never });

function createList(fields: Partial<ParticleList> & Pick<ParticleList, 'count' | 'targets' | 'coordinates'>) {
  return {
    targetInfo: new Map(),
    getParticleLabel: (i: number) => `particle ${i + 1}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    ...fields,
  } as ParticleList;
}

/** three particles of two targets (0, 1), entities on two of them, rotated alike but the last */
function createTemplate() {
  const r = Quat.setAxisAngle(Quat(), Vec3.create(1, 0, 0), Math.PI / 3);
  return createList({
    label: 'template',
    count: 3,
    coordinates: Float32Array.of(1, 0, 0, 0, 2, 0, 0, 0, 3),
    rotations: Float32Array.of(...r, ...r, 0, 0, 0, 1),
    targets: Int32Array.of(0, 0, 1),
    targetInfo: new Map([
      [0, { name: 'a' }],
      [1, { name: 'b' }],
    ]),
    targetMapping: new Map([
      [0, leaf('a')],
      [1, leaf('b')],
    ]),
    entities: Int32Array.of(0, -1, 1),
    entityInfo: new Map([
      [0, { name: 'A' }],
      [1, { name: 'B' }],
    ]),
    radii: Float32Array.of(0.5, 0.5, 1),
    getParticleLabel: (i) => `template ${i + 1}`,
  });
}

/** four particles: three placing the template (target 2, centered at 0), turned about z or not, one a leaf (target 0) */
/** `center` null: none given (the template's boundary center) */
function createTop(template: ParticleList, center: Vec3 | null = Vec3.create(0, 0, 0)) {
  const turn = Quat.setAxisAngle(Quat(), Vec3.create(0, 0, 1), Math.PI / 2);
  const id = Quat.identity();
  return createList({
    label: 'top',
    count: 4,
    coordinates: Float32Array.of(10, 0, 0, 0, 10, 0, 5, 5, 5, 0, 0, 20),
    rotations: Float32Array.of(...id, ...turn, ...id, ...id),
    targets: Int32Array.of(2, 2, 0, 2),
    targetInfo: new Map([
      [0, { name: 'leaf' }],
      [2, { name: 'copy' }],
    ]),
    targetMapping: new Map<number, ParticleTarget>([
      [0, leaf('top')],
      [2, { kind: 'particles', particles: template, center: center ?? undefined }],
    ]),
    entities: Int32Array.of(0, 0, 0, 0),
    entityInfo: new Map([[0, { name: 'T' }]]),
    radii: Float32Array.of(5, 5, 2, 5),
  });
}

const position = (list: ParticleList, i: number) => GridPositions.get(list.coordinates, i, Vec3());

describe('particle nesting', () => {
  it('expands a template at each particle placing it: paths and indices round-trip', () => {
    const top = createTop(createTemplate());
    const { count, offsets } = ParticleNesting.getExpansion(top);
    expect(count).toBe(3 + 3 + 1 + 3);
    expect(Array.from(offsets!)).toEqual([0, 3, 6, 7, 10]);
    const path: ParticlePath = { top: 0, local: 0 };
    const paths: [number, number][] = [];
    for (let i = 0; i < count; ++i) {
      ParticleNesting.getPath(top, i, path);
      paths.push([path.top, path.local]);
      expect(ParticleNesting.getExpandedIndex(top, path.top, path.local)).toBe(i);
    }
    expect(paths).toEqual([
      [0, 0],
      [0, 1],
      [0, 2],
      [1, 0],
      [1, 1],
      [1, 2],
      [2, 0],
      [3, 0],
      [3, 1],
      [3, 2],
    ]);
  });

  it('flattens: poses composed, ids moved past the top list’s, entities inherited where none, labels by path', () => {
    const template = createTemplate();
    const top = createTop(template);
    const flat = ParticleNesting.flatten(top);
    expect(ParticleNesting.flatten(top)).toBe(flat);
    expect(flat.count).toBe(10);
    // the second copy turned 90° about z: (1, 0, 0) to (0, 1, 0), at (0, 10, 0)
    expect(Array.from(position(flat, 3))).toEqual([expect.closeTo(0, 5), expect.closeTo(11, 5), expect.closeTo(0, 5)]);
    expect(Array.from(position(flat, 4))).toEqual([expect.closeTo(-2, 5), expect.closeTo(10, 5), expect.closeTo(0, 5)]);
    // every instance as its top pose composed with its template particle's
    const q = Quat();
    const qc = Quat();
    const p = Vec3();
    for (let i = 0; i < flat.count; ++i) {
      const { top: t, local } = ParticleNesting.getPath(top, i, { top: 0, local: 0 });
      const rotation = ParticleRotations.get(flat.rotations, i, Quat());
      if (t === 2) {
        expect(Array.from(position(flat, i))).toEqual([5, 5, 5]);
        continue;
      }
      ParticleRotations.get(top.rotations, t, q);
      Vec3.add(p, Vec3.transformQuat(p, position(template, local), q), position(top, t));
      expect(Vec3.distance(position(flat, i), p)).toBeLessThan(1e-5);
      Quat.multiply(qc, q, ParticleRotations.get(template.rotations, local, qc));
      // q and -q the same turn
      expect(Math.abs(Quat.dot(rotation, qc))).toBeCloseTo(1, 6);
    }
    // the leaf keeps its target; the template's targets past the top list's ids (0 and 2 taken): 3 and 4
    const targets = Array.from({ length: 10 }, (_, i) => ParticleValues.get(flat.targets, i));
    expect(targets).toEqual([3, 3, 4, 3, 3, 4, 0, 3, 3, 4]);
    expect(flat.targetMapping!.get(0)).toBe(top.targetMapping!.get(0));
    expect(flat.targetMapping!.get(3)).toBe(template.targetMapping!.get(0));
    expect(flat.targetMapping!.get(4)).toBe(template.targetMapping!.get(1));
    expect(flat.targetMapping!.has(2)).toBe(false);
    expect([...flat.targetInfo.entries()]).toEqual([
      [0, { name: 'leaf' }],
      [3, { name: 'a' }],
      [4, { name: 'b' }],
    ]);
    // the template's entities past the top's (0): 1, 2; its particle without one takes its top particle's (0)
    const entities = Array.from({ length: 10 }, (_, i) => ParticleValues.get(flat.entities!, i));
    expect(entities).toEqual([1, 0, 2, 1, 0, 2, 0, 1, 0, 2]);
    expect(flat.entityInfo!.get(2)).toEqual({ name: 'B' });
    expect(Array.from({ length: 10 }, (_, i) => ParticleValues.get(flat.radii!, i))).toEqual([
      0.5, 0.5, 1, 0.5, 0.5, 1, 2, 0.5, 0.5, 1,
    ]);
    expect(flat.getParticleLabel(4)).toBe('particle 2 › template 2');
    expect(flat.getParticleLabel(6)).toBe('particle 3');
  });

  it('places a template about its center: given, or its boundary’s by default', () => {
    const template = createTemplate();
    const center = Particle.getBoundary(template).sphere.center;
    const flat = ParticleNesting.flatten(createTop(template, null));
    // the first copy unturned at (10, 0, 0): a particle at its position less the center, moved there
    const expected = Vec3.add(Vec3(), Vec3.sub(Vec3(), position(template, 0), center), Vec3.create(10, 0, 0));
    expect(Vec3.distance(position(flat, 0), expected)).toBeLessThan(1e-5);
  });

  it('bounds a list with templates by its copies: every instance inside, as tight as the flattened list’s', () => {
    const top = createTop(createTemplate());
    const flat = ParticleNesting.flatten(top);
    const { sphere, box } = Particle.getBoundary(top);
    for (let i = 0; i < flat.count; ++i) {
      const p = position(flat, i);
      const r = ParticleValues.get(flat.radii!, i);
      expect(Vec3.distance(p, sphere.center) + r).toBeLessThanOrEqual(sphere.radius + 1e-4);
      for (let k = 0; k < 3; ++k) {
        expect(p[k] - r).toBeGreaterThanOrEqual(box.min[k] - 1e-4);
        expect(p[k] + r).toBeLessThanOrEqual(box.max[k] + 1e-4);
      }
    }
    // a copy's sphere holds its template's, so at most a little larger than the instances' own
    expect(sphere.radius).toBeLessThan(Particle.getBoundary(flat).sphere.radius * 1.25);
  });

  it('refuses nesting deeper than a level, a list placing itself, a slot list as a template', () => {
    const template = createTemplate();
    const nested = createTop(template);
    ParticleNesting.check(nested);
    expect(() => ParticleNesting.check(createTop(nested))).toThrow(/templates/);
    const self = createTop(template) as { targetMapping: Map<number, ParticleTarget> };
    self.targetMapping.set(2, { kind: 'particles', particles: self as unknown as ParticleList });
    expect(() => ParticleNesting.check(self as unknown as ParticleList)).toThrow(/own list/);
    const paged = { ...template, pages: { pages: [] } } as unknown as ParticleList;
    expect(() => ParticleNesting.check(createTop(paged))).toThrow(/slot list/);
  });

  it('leaves a list without templates as it is: an instance a particle', () => {
    const template = createTemplate();
    expect(ParticleNesting.hasTemplates(template)).toBe(false);
    expect(ParticleNesting.flatten(template)).toBe(template);
    expect(ParticleNesting.getExpansion(template)).toEqual({ count: 3 });
    expect(ParticleNesting.getPath(template, 2, { top: 0, local: 0 })).toEqual({ top: 2, local: 0 });
    expect(ParticleNesting.getExpandedIndex(template, 2, 0)).toBe(2);
  });

  /** deterministic subsets of the flattened list's instances (a few whole copies among them) */
  function subsets(count: number, n: number) {
    let seed = 7;
    const random = () => ((seed = (seed * 16807) % 2147483647) - 1) / 2147483646;
    const out: number[][] = [[], Array.from({ length: count }, (_, i) => i)];
    for (let k = 0; k < n; ++k) {
      const p = random();
      out.push(Array.from({ length: count }, (_, i) => i).filter(() => random() < p));
    }
    return out;
  }

  it('holds instances by path: a copy, an instance, a template particle in every copy; sizes, flat indices', () => {
    const template = createTemplate();
    const top = createTop(template);
    const flatIndices = (loci: Particle.Loci) => OrderedSet.toArray(ParticleNesting.toFlatLoci(loci).indices);
    // copy 1 (its 3 instances), the leaf (one), instance 2 of copy 3
    expect(flatIndices(ParticleNesting.ofCopy(top, 1))).toEqual([3, 4, 5]);
    expect(flatIndices(ParticleNesting.ofCopy(top, 2))).toEqual([6]);
    expect(flatIndices(ParticleNesting.ofInstance(top, 3, 2))).toEqual([9]);
    expect(Particle.lociSize(ParticleNesting.ofCopy(top, 1))).toBe(3);
    // template particle 0 in every copy: [*, 0]
    const all0 = ParticleNesting.ofTemplateParticles(top, template, OrderedSet.ofSingleton(0));
    expect(flatIndices(all0)).toEqual([0, 3, 7]);
    expect(Array.from(all0.partial!.tops)).toEqual([0, 1, 3]);
    // all of a copy's template particles: the copy whole, no partial
    const allOf1 = ParticleNesting.createLoci(top, OrderedSet.Empty, [1], [OrderedSet.ofRange(0, 2)]);
    expect(Particle.areLociEqual(allOf1, ParticleNesting.ofCopy(top, 1))).toBe(true);
    expect(allOf1.partial).toBeUndefined();
    // the whole list: an interval of tops, its flat loci the flattened list's every instance
    const whole = Particle.Loci(top, OrderedSet.ofBounds(0, 4));
    expect(ParticleNesting.toFlatLoci(whole).indices).toEqual(OrderedSet.ofBounds(0, 10));
    expect(ParticleNesting.toFlatLoci(whole).particles).toBe(ParticleNesting.flatten(top));
  });

  it('unites, intersects and subtracts by path as the flattened list’s indices do; from flat indices and back', () => {
    const top = createTop(createTemplate());
    const sets = subsets(10, 24).map((xs) => ParticleNesting.fromFlatIndices(top, OrderedSet.ofSortedArray(xs)));
    const flat = (loci: Particle.Loci) => ParticleNesting.toFlatLoci(loci).indices;
    for (const a of sets) {
      // canonical: from its own flat indices the same loci
      expect(Particle.areLociEqual(ParticleNesting.fromFlatIndices(top, flat(a)), a)).toBe(true);
      expect(Particle.lociSize(a)).toBe(OrderedSet.size(flat(a)));
      for (const b of sets.slice(0, 8)) {
        const ops = [
          [ParticleNesting.union(a, b), OrderedSet.union(flat(a), flat(b))],
          [ParticleNesting.intersect(a, b), OrderedSet.intersect(flat(a), flat(b))],
          [ParticleNesting.subtract(a, b), OrderedSet.subtract(flat(a), flat(b))],
        ] as const;
        for (const [loci, expected] of ops) {
          expect(OrderedSet.toArray(flat(loci))).toEqual(OrderedSet.toArray(expected));
          expect(Particle.areLociEqual(loci, ParticleNesting.fromFlatIndices(top, expected))).toBe(true);
        }
      }
    }
    expect(() => ParticleNesting.union(sets[2], Particle.Loci(createTop(createTemplate()), OrderedSet.Empty))).toThrow(
      /different lists/,
    );
  });

  it('queries by level: the top particles, a template’s particles in every copy, both intersected', () => {
    const template = createTemplate();
    const top = createTop(template);
    const flat = ParticleNesting.flatten(top);
    // the template's entity B (its id 1) in every copy: the flattened list's entity 2 (moved past the top's)
    const byEntity = ParticleNesting.selectInTemplate(
      top,
      template,
      (j) => ParticleValues.get(template.entities!, j) === 1,
    );
    const expected = Array.from({ length: flat.count }, (_, i) => i).filter(
      (i) => ParticleValues.get(flat.entities!, i) === 2,
    );
    expect(OrderedSet.toArray(ParticleNesting.toFlatLoci(byEntity).indices)).toEqual(expected);
    // copies turned (a top's rotation not the identity): copy 1 alone
    const turned = ParticleNesting.selectTop(top, (i) => ParticleRotations.get(top.rotations, i, Quat())[3] < 0.99);
    expect(OrderedSet.toArray(turned.indices)).toEqual([1]);
    const both = ParticleNesting.intersect(byEntity, turned);
    expect(OrderedSet.toArray(ParticleNesting.toFlatLoci(both).indices)).toEqual([5]);
  });

  it('bounds, labels and remaps a loci by path', () => {
    const template = createTemplate();
    const top = createTop(template);
    const flat = ParticleNesting.flatten(top);
    for (const xs of subsets(10, 12).slice(2)) {
      if (xs.length === 0) continue;
      const loci = ParticleNesting.fromFlatIndices(top, OrderedSet.ofSortedArray(xs));
      const sphere = Particle.getBoundingSphere(loci);
      for (const i of xs) {
        const r = ParticleValues.get(flat.radii!, i);
        expect(Vec3.distance(position(flat, i), sphere.center) + r).toBeLessThanOrEqual(sphere.radius + 1e-4);
      }
    }
    expect(Particle.getLabel(ParticleNesting.ofInstance(top, 1, 1))).toBe('particle 2 › template 2');
    expect(Particle.getLabel(ParticleNesting.ofCopy(top, 1))).toBe('particle 2 · 3 particles');
    expect(Particle.getLabel(ParticleNesting.ofCopy(top, 2))).toBe('particle 3');
    expect(Particle.getLabel(Particle.Loci(top, OrderedSet.ofBounds(0, 4)))).toBe('10 particles');
    expect(Particle.isLociEmpty(ParticleNesting.createLoci(top, OrderedSet.Empty, [0], [OrderedSet.Empty]))).toBe(true);
    // to a list of the first two top particles: the others' parts dropped
    const shorter = { ...top, count: 2, _propertyData: Object.create(null) } as ParticleList;
    const loci = ParticleNesting.union(ParticleNesting.ofInstance(top, 0, 1), ParticleNesting.ofInstance(top, 3, 0));
    const remapped = Particle.remapLoci(loci, shorter);
    expect(Array.from(remapped.partial!.tops)).toEqual([0]);
  });

  it('bounds a template particle set in every copy by its sphere once, placed: every instance inside', () => {
    // 40 template particles on a helix off the copies' turn axis, every other one held in each copy
    const n = 40;
    const coordinates = new Float32Array(n * 3);
    for (let i = 0; i < n; ++i) coordinates.set([10 + Math.cos(i * 0.5) * 4, Math.sin(i * 0.5) * 4, i * 0.3], i * 3);
    const template = createList({
      count: n,
      coordinates,
      targets: new Int32Array(n),
      radii: new Float32Array(n).fill(1),
    });
    const top = createTop(template);
    const flat = ParticleNesting.flatten(top);
    const half = OrderedSet.ofSortedArray(Int32Array.from({ length: n / 2 }, (_, k) => k * 2));
    const loci = ParticleNesting.ofTemplateParticles(top, template, half);
    const sphere = Particle.getBoundingSphere(loci);
    let reach = 0;
    ParticleNesting.forEachInstance(loci, (t, local) => {
      const i = ParticleNesting.getExpandedIndex(top, t, local);
      const d = Vec3.distance(position(flat, i), sphere.center) + 1;
      expect(d).toBeLessThanOrEqual(sphere.radius + 1e-4);
      reach = Math.max(reach, d);
    });
    // a sphere of spheres: looser than the instances', not by much
    expect(sphere.radius).toBeLessThan(reach * 1.25);
  });

  it('gives a template its leaf list: its particles in the flattened ids, each copy’s instances their values', () => {
    const template = createTemplate();
    const top = createTop(template);
    const flat = ParticleNesting.flatten(top);
    const ids = ParticleNesting.getIdSpace(top);
    expect(ids.targetInfo).toBe(flat.targetInfo);
    const leaves = ParticleNesting.getLeafList(top, template);
    expect(ParticleNesting.getLeafList(top, template)).toBe(leaves);
    expect(leaves.count).toBe(template.count);
    for (const t of ParticleNesting.getCopiesOf(top, template)) {
      for (let j = 0; j < template.count; ++j) {
        const i = ParticleNesting.getExpandedIndex(top, t, j);
        expect(ParticleValues.get(leaves.targets, j)).toBe(ParticleValues.get(flat.targets, i));
        // an entity of its own the flattened one's; none its own (-1), the flattened list's the copy's top particle's
        const e = ParticleValues.get(leaves.entities!, j);
        if (e >= 0) expect(e).toBe(ParticleValues.get(flat.entities!, i));
        else expect(ParticleValues.get(flat.entities!, i)).toBe(ParticleValues.get(top.entities!, t));
      }
    }
    expect(leaves.targetMapping!.get(ParticleValues.get(leaves.targets, 2))).toBe(template.targetMapping!.get(1));
    expect(() => ParticleNesting.getLeafList(top, createTemplate())).toThrow(/not a template/);
  });

  it('places its templates’ density grids by their copies: inside a copy as its template’s, none on its edges', () => {
    const template = createTemplate();
    // a 7 Å cube of 7³ blocks, each of density 1
    ParticleDensityGrid.set(template, {
      min: Vec3.create(0, 0, 0),
      max: Vec3.create(7, 7, 7),
      dimensions: [7, 7, 7],
      density: new Float32Array(343).fill(1),
    });
    const top = createTop(template);
    // the list's blocks about 2 Å (its 27 Å along z in 13), coarser than the template's (a list of many copies' are),
    // the copies' edges within blocks (7 Å not a multiple of them)
    const grid = ParticleNesting.getDensityGrid(top, 13)!;
    expect(ParticleNesting.getDensityGrid(top, 13)).toBe(grid);
    const [bx, by, bz] = grid.dimensions;
    const s = [0, 1, 2].map((a) => (grid.max[a] - grid.min[a]) / grid.dimensions[a]);
    // the copies' cubes (about the origin): at 10, 0, 0; turned 90° about z at 0, 10, 0; at 0, 0, 20
    const cubes = [
      [10, 17, 0, 7, 0, 7],
      [-7, 0, 10, 17, 0, 7],
      [0, 7, 0, 7, 20, 27],
    ];
    const within = (x: number, y: number, z: number) =>
      cubes.some(
        (c) =>
          x >= c[0] - 1e-6 &&
          x <= c[1] + 1e-6 &&
          y >= c[2] - 1e-6 &&
          y <= c[3] + 1e-6 &&
          z >= c[4] - 1e-6 &&
          z <= c[5] + 1e-6,
      );
    let total = 0;
    for (let z = 0, b = 0; z < bz; ++z) {
      for (let y = 0; y < by; ++y) {
        for (let x = 0; x < bx; ++x, ++b) {
          const v = grid.density[b];
          if (!(v > 0)) continue;
          total += v * s[0] * s[1] * s[2];
          // a block with density lies whole within a copy (its corners): none on a copy's edge, a gap open
          const x0 = grid.min[0] + x * s[0];
          const y0 = grid.min[1] + y * s[1];
          const z0 = grid.min[2] + z * s[2];
          expect(within(x0, y0, z0) && within(x0 + s[0], y0 + s[1], z0 + s[2])).toBe(true);
        }
      }
    }
    // each copy's inside kept (the blocks on its edges dropped): some of it, less than the whole
    expect(total).toBeGreaterThan(3 * 3 * 3 * 3);
    expect(total).toBeLessThan(3 * 343);
    // a copy's center its template's density, between the copies none
    const at = (p: Vec3) => {
      const i = Math.floor((p[0] - grid.min[0]) / s[0]);
      const j = Math.floor((p[1] - grid.min[1]) / s[1]);
      const k = Math.floor((p[2] - grid.min[2]) / s[2]);
      return grid.density[i + bx * (j + by * k)];
    };
    expect(at(Vec3.create(13.5, 3.5, 3.5))).toBeGreaterThan(0.5);
    expect(at(Vec3.create(-3.5, 13.5, 3.5))).toBeGreaterThan(0.5);
    expect(at(Vec3.create(5, 5, 12))).toBe(0);
  });

  it('samples a list by target: every 8^L-th of a target, at least one, spread over its box, its area kept', () => {
    // a 10³ grid of target 0 and seven particles of target 1 in a corner
    const n = 1000 + 7;
    const coordinates = new Float32Array(n * 3);
    for (let i = 0; i < 1000; ++i) coordinates.set([i % 10, Math.floor(i / 10) % 10, Math.floor(i / 100)], i * 3);
    for (let i = 0; i < 7; ++i) coordinates.set([0.5, 0.5, i * 0.1], (1000 + i) * 3);
    const list = createList({
      count: n,
      coordinates,
      targets: Int32Array.from({ length: n }, (_, i) => (i < 1000 ? 0 : 1)),
    });
    const sample = Particle.sampleByTarget(list, 1)!;
    const { list: s, indices } = sample;
    expect(s.count).toBe(125 + 1);
    expect(new Set(indices).size).toBe(s.count);
    // a run a target, its scale its area's: kept times scale squared its particles
    const scales = s.scales!;
    expect(ParticleValues.isRuns(scales)).toBe(true);
    let kept0 = 0;
    for (let i = 0; i < s.count; ++i) {
      const t = ParticleValues.get(s.targets, i);
      expect(t).toBe(ParticleValues.get(list.targets, indices[i]));
      expect(position(s, i)).toEqual(position(list, indices[i]));
      if (t === 0) kept0 += 1;
      const seen = t === 0 ? 1000 : 7;
      const kept = t === 0 ? 125 : 1;
      expect(ParticleValues.get(scales, i)).toBeCloseTo(Math.sqrt(seen / kept), 5);
    }
    expect(kept0).toBe(125);
    // spread over the grid: its sample in every octant of it
    const octants = new Set<number>();
    for (let i = 0; i < s.count; ++i) {
      if (ParticleValues.get(s.targets, i) !== 0) continue;
      const p = position(s, i);
      octants.add((p[0] >= 5 ? 1 : 0) + (p[1] >= 5 ? 2 : 0) + (p[2] >= 5 ? 4 : 0));
    }
    expect(octants.size).toBe(8);
    // deterministic; a list of scales by run of its own none
    expect(Particle.sampleByTarget(list, 1)!.indices).toEqual(indices);
    expect(Particle.sampleByTarget(s, 1)).toBeUndefined();
    // a level from the one before (every eighth of a target's kept) as from all the particles
    const order = Particle.getSampleOrder(list);
    const second = Particle.sampleByTarget(list, 2, order)!;
    const fromFirst = Particle.sampleByTarget(list, 2, order, sample)!;
    expect(fromFirst.indices).toEqual(second.indices);
    expect(fromFirst.list.scales).toEqual(second.list.scales);
    expect(second.list.count).toBe(16 + 1);
  });

  it("samples a list stored in sample order by its own order: a target's every 8^L-th in its run, or as stored", () => {
    // 20 particles of target 0 then 10 of target 1 (runs), stored in sample order (no sort)
    const n = 30;
    const runs = (offsets: number[], values: number[]) =>
      ParticleValues.ofRuns(Int32Array.from(offsets), Int32Array.from(values));
    const ordered = (targets: ParticleList['targets']) => {
      const list = createList({ count: n, coordinates: new Float32Array(n * 3).map((_, k) => k), targets });
      Particle.setSampleOrdered(list);
      return list;
    };
    const a = ordered(runs([0, 20, 30], [0, 1]));
    expect(Array.from(Particle.getSampleOrder(a))).toEqual(Array.from({ length: n }, (_, i) => i));
    const s = Particle.sampleByTarget(a, 1)!;
    expect(Array.from(s.indices)).toEqual([0, 8, 16, 20, 28]);
    expect(ParticleValues.get(s.list.scales!, 0)).toBeCloseTo(Math.sqrt(20 / 3), 5);
    // a target in two runs apart: picked as the order gives it (every 8th of its particles as stored, all its runs)
    const b = ordered(runs([0, 10, 20, 30], [0, 1, 0]));
    expect(Array.from(Particle.sampleByTarget(b, 1)!.indices)).toEqual([0, 8, 26, 10, 18]);
  });

  it("gives a sample its targets' entities by run where a target's particles have one, else a particle's", () => {
    const n = 64;
    const coordinates = new Float32Array(n * 3).map((_, k) => k);
    const targets = Int32Array.from({ length: n }, (_, i) => (i < 32 ? 0 : 1));
    // target 0's particles of one entity, target 1's of two (alternating)
    const entities = Int32Array.from({ length: n }, (_, i) => (i < 32 ? 5 : 6 + (i % 2)));
    const list = createList({ count: n, coordinates, targets, entities });
    const s = Particle.sampleByTarget(list, 1)!;
    for (let k = 0; k < s.list.count; ++k) {
      expect(ParticleValues.get(s.list.entities!, k)).toBe(entities[s.indices[k]]);
      expect(ParticleValues.get(s.list.targets, k)).toBe(targets[s.indices[k]]);
    }
    // a list of one entity a target: its sample's entities runs (a value a target)
    const uniform = createList({ count: n, coordinates, targets, entities: Int32Array.from(targets, (t) => 9 + t) });
    const u = Particle.sampleByTarget(uniform, 1)!;
    expect(ParticleValues.isRuns(u.list.entities!)).toBe(true);
    for (let k = 0; k < u.list.count; ++k)
      expect(ParticleValues.get(u.list.entities!, k)).toBe(9 + targets[u.indices[k]]);
  });

  it('gives a template its coarse levels: sampled while a level keeps half the one before, eight a target', () => {
    // a template of 20³ particles of two targets, placed twice
    const n = 8000;
    const coordinates = new Float32Array(n * 3);
    for (let i = 0; i < n; ++i) coordinates.set([i % 20, Math.floor(i / 20) % 20, Math.floor(i / 400)], i * 3);
    const big = createList({
      label: 'big',
      count: n,
      coordinates,
      targets: Int32Array.from({ length: n }, (_, i) => i % 2),
      targetInfo: new Map([
        [0, { name: 'a' }],
        [1, { name: 'b' }],
      ]),
      targetMapping: new Map([
        [0, leaf('a')],
        [1, leaf('b')],
      ]),
    });
    const top = createTop(big);
    const levels = ParticleNesting.getTemplateLevels(top, big);
    // the fourth's one a target not kept (eight a target on average at least)
    expect(levels.map((l) => l.level)).toEqual([1, 2, 3]);
    expect(levels.map((l) => l.list.count)).toEqual([1000, 126, 16]);
    // the leaf list's ids (the flattened list's), cached
    const leaves = ParticleNesting.getLeafList(top, big);
    for (const { list, indices } of levels) {
      for (let i = 0; i < list.count; ++i) {
        expect(ParticleValues.get(list.targets, i)).toBe(ParticleValues.get(leaves.targets, indices[i]));
      }
    }
    expect(ParticleNesting.getTemplateLevels(top, big)).toBe(levels);
    // a template of a few particles a target: none
    const small = createTemplate();
    expect(ParticleNesting.getTemplateLevels(createTop(small), small)).toEqual([]);
  });

  it("gives a list other targets with its targets' caches made again, its particles' kept", () => {
    const template = createTemplate();
    const top = createTop(template);
    // cached of its targets: its expansion, id space, flattened list, boundary; of its particles: its transforms
    expect(ParticleNesting.getExpandedCount(top)).toBe(10);
    const flat = ParticleNesting.flatten(top);
    const transforms = getParticleTransforms(top);
    Particle.getBoundary(top);
    // another template at the copies' target: four particles, not three
    const other = createList({
      ...createTemplate(),
      count: 4,
      coordinates: new Float32Array(12),
      targets: new Int32Array(4),
    });
    const mapping = new Map(top.targetMapping!);
    mapping.set(2, { kind: 'particles', particles: other, center: Vec3.create(0, 0, 0) });
    const retargeted = Particle.withTargets(top, mapping);
    expect(ParticleNesting.getExpandedCount(retargeted)).toBe(4 + 4 + 1 + 4);
    expect(ParticleNesting.flatten(retargeted)).not.toBe(flat);
    expect([...ParticleNesting.getIdSpace(retargeted).moves.keys()]).toEqual([other]);
    expect(getParticleTransforms(retargeted)).toBe(transforms);
    // a flat list's boundary kept (its particles' alone)
    const leaves = createTemplate();
    const boundary = Particle.getBoundary(leaves);
    expect(Particle.getBoundary(Particle.withTargets(leaves, new Map()))).toBe(boundary);
  });
});
