/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { CustomProperties } from '../../custom-property.js';
import { Particle, type ParticleList, ParticleValues } from '../particle-list.js';
import { type ParticlePage, ParticleSlots } from '../particle-slots.js';

const Count = 30;

function createList(): ParticleList {
  const coordinates = Float32Array.from({ length: Count * 3 }, (_, k) => k);
  return {
    label: 'list',
    count: Count,
    targets: Int32Array.from({ length: Count }, (_, i) => (i < 15 ? 1 : 2)),
    targetInfo: new Map([
      [1, { name: 'one' }],
      [2, { name: 'two' }],
    ]),
    entities: Int32Array.from({ length: Count }, (_, i) => i % 3),
    coordinates,
    getParticleLabel: (i) => `particle ${i}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

const list = createList();
/** `[a, b)` */
const range = (a: number, b: number) => Array.from({ length: b - a }, (_, i) => a + i);
/** particles `[a, b)` of the list as a page at slot `start` */
const pageOf = (start: number, a: number, b: number): ParticlePage => ({
  start,
  list: Particle.subset(
    list,
    Int32Array.from({ length: b - a }, (_, i) => a + i),
  ),
});

describe('particle slots', () => {
  it('reads a flat list as one page of all its particles', () => {
    const pages: ParticlePage[] = [];
    ParticleSlots.forEachPage(list, (p) => pages.push(p));
    expect(pages).toEqual([{ start: 0, list }]);
    expect(ParticleSlots.isPaged(list)).toBe(false);
    expect(ParticleSlots.getLiveCount(list)).toBe(Count);
    expect(ParticleSlots.findPage(list, 29)).toBe(0);
    expect(ParticleSlots.findPage(list, 30)).toBe(-1);
    expect(ParticleSlots.locate(list, 7)).toEqual({ list, index: 7 });
    expect(ParticleSlots.flatten(list)).toBe(list);
  });

  it('holds its particles in pages of slots, the slots between them gaps', () => {
    const a = pageOf(0, 0, 10);
    const b = pageOf(100, 10, 15);
    const c = pageOf(40, 15, 30);
    const slots = ParticleSlots.create(list, 200, [a, b, c]);
    expect(slots.count).toBe(200);
    expect(ParticleSlots.isPaged(slots)).toBe(true);
    expect(ParticleSlots.getLiveCount(slots)).toBe(30);
    // by start
    expect(slots.pages!.pages).toEqual([a, c, b]);
    expect([0, 9, 10, 39, 40, 54, 55, 100, 104, 105].map((s) => ParticleSlots.findPage(slots, s))).toEqual([
      0, 0, -1, -1, 1, 1, -1, 2, 2, -1,
    ]);
    const s = ParticleSlots.locate(slots, 42)!;
    expect(s.list).toBe(c.list);
    expect(s.index).toBe(2);
    expect(ParticleSlots.locate(slots, 50 + 10)).toBeUndefined();
    // a slot's label its particle's, a gap's its own
    expect(slots.getParticleLabel(42)).toBe('particle 17');
    expect(slots.getParticleLabel(70)).toBe('Slot 71 (empty)');
    // what the pages share the list's own
    expect(slots.targetInfo).toBe(list.targetInfo);
    // a first version: no changes
    expect(slots.pagesFrom).toBeUndefined();
  });

  it('flattens to a flat list by slot, a particle keyed by its slot', () => {
    const slots = ParticleSlots.create(list, 200, [pageOf(100, 10, 15), pageOf(0, 0, 10), pageOf(40, 15, 30)]);
    const flat = ParticleSlots.flatten(slots);
    expect(flat.count).toBe(30);
    expect(flat.pages).toBeUndefined();
    const slotOf = (i: number) => flat.keys![i];
    expect(Array.from({ length: 30 }, (_, i) => slotOf(i))).toEqual([
      ...Array.from({ length: 10 }, (_, i) => i),
      ...Array.from({ length: 15 }, (_, i) => 40 + i),
      ...Array.from({ length: 5 }, (_, i) => 100 + i),
    ]);
    const p = Vec3();
    const q = Vec3();
    for (let i = 0; i < flat.count; ++i) {
      const s = ParticleSlots.locate(slots, slotOf(i))!;
      expect(Array.from(GridPositions.get(flat.coordinates, i, p))).toEqual(
        Array.from(GridPositions.get(s.list.coordinates, s.index, q)),
      );
      expect(ParticleValues.get(flat.targets, i)).toBe(ParticleValues.get(s.list.targets, s.index));
      expect(ParticleValues.get(flat.entities!, i)).toBe(ParticleValues.get(s.list.entities!, s.index));
      expect(flat.getParticleLabel(i)).toBe(s.list.getParticleLabel(s.index));
    }
  });

  it('adds and removes pages, the others and their slots kept, the changes said', () => {
    const a = pageOf(0, 0, 10);
    const b = pageOf(40, 10, 20);
    const v1 = ParticleSlots.create(list, 100, [a, b]);
    const c = pageOf(80, 20, 30);
    const v2 = ParticleSlots.update(v1, [c], [0]);
    expect(v2.pages!.pages).toEqual([b, c]);
    expect(v2.pages!.pages[0]).toBe(b);
    expect(ParticleSlots.locate(v2, 45)!.list).toBe(b.list);
    expect(ParticleSlots.locate(v2, 5)).toBeUndefined();
    const changes = ParticleSlots.getChanges(v1, v2)!;
    expect(changes.added).toEqual([c]);
    expect(changes.removed).toEqual([a]);
    expect(changes.posed).toEqual([]);
    // a page replaced by its particles posed otherwise: posed
    const moved = { start: 40, list: Particle.withPoses(b.list, new Float32Array(b.list.count * 3)) };
    const posed = ParticleSlots.getChanges(v2, ParticleSlots.update(v2, [moved], [40]))!;
    expect(posed.posed).toEqual([moved]);
    expect(posed.added).toEqual([]);
    expect(posed.removed).toEqual([]);
    // one step only: from the first to a third, built again
    const v3 = ParticleSlots.update(v2, [a]);
    expect(ParticleSlots.getChanges(v2, v3)!.added).toEqual([a]);
    expect(ParticleSlots.getChanges(v1, v3)).toBeUndefined();
    expect(ParticleSlots.getLiveCount(v3)).toBe(30);
  });

  it('holds no earlier version of itself: a version reaches one step back, not every version and page since the first', () => {
    const a = pageOf(0, 0, 10);
    const v1 = ParticleSlots.create(list, 100, [a]);
    const v2 = ParticleSlots.update(v1, [pageOf(40, 10, 20)], [0]);
    const v3 = ParticleSlots.update(v2, [pageOf(80, 20, 30)], [40]);
    // what v3 references, walked: neither list before it, nor page a (removed by v2)
    const seen = new Set<unknown>();
    const walk = (o: unknown) => {
      if (o === null || typeof o !== 'object' || seen.has(o) || ArrayBuffer.isView(o)) return;
      seen.add(o);
      for (const v of Object.values(o)) walk(v);
    };
    walk(v3.pagesFrom);
    expect(seen.has(v1)).toBe(false);
    expect(seen.has(v2)).toBe(false);
    expect(seen.has(a)).toBe(false);
    expect(ParticleSlots.getChanges(v2, v3)!.removed.map((p) => p.start)).toEqual([40]);
  });

  it("keeps the sphere all its slots' particles are within over its versions, whatever pages are live", () => {
    const sphere = Sphere3D.create(Vec3.create(1, 2, 3), 50);
    const v1 = ParticleSlots.create(list, 100, [pageOf(0, 0, 10)], undefined, sphere);
    expect(v1.pages!.sphere).toBe(sphere);
    const v2 = ParticleSlots.update(v1, [pageOf(40, 10, 20)], [0]);
    expect(v2.pages!.sphere).toBe(sphere);
    expect(ParticleSlots.create(list, 100, []).pages!.sphere).toBeUndefined();
  });

  it('refuses pages overlapping, past the slots or of pages', () => {
    expect(() => ParticleSlots.create(list, 100, [pageOf(0, 0, 10), pageOf(5, 10, 20)])).toThrow(/overlapping/);
    expect(() => ParticleSlots.create(list, 15, [pageOf(10, 0, 10)])).toThrow(/past the slots/);
    const slots = ParticleSlots.create(list, 100, [pageOf(0, 0, 10)]);
    expect(() => ParticleSlots.create(list, 100, [{ start: 0, list: slots }])).toThrow(/of pages/);
    expect(() => ParticleSlots.update(list, [])).toThrow(/not a slot list/);
  });
});

describe('particle slots: locations, loci and bounds', () => {
  const pages = () => [pageOf(0, 0, 10), pageOf(40, 10, 20), pageOf(80, 20, 30)];

  it("places a slot's location by its page's particle, a gap at the origin", () => {
    const slots = ParticleSlots.create(list, 100, pages());
    const p = Vec3();
    const q = Vec3();
    Particle.position(p, Particle.Location(slots, 45));
    GridPositions.get(list.coordinates, 15, q);
    expect(Array.from(p)).toEqual(Array.from(q));
    expect(Array.from(Particle.position(p, Particle.Location(slots, 60)))).toEqual([0, 0, 0]);
  });

  it("bounds a slot list by its pages' bounds and a loci of slots by their particles", () => {
    const slots = ParticleSlots.create(list, 100, pages());
    const whole = Particle.getBoundary(list).sphere;
    const b = Particle.getBoundary(slots).sphere;
    // its pages' spheres' union: around the whole's
    expect(b.radius).toBeGreaterThanOrEqual(whole.radius - 1e-3);
    expect(Vec3.distance(b.center, whole.center)).toBeLessThan(whole.radius);
    // a loci of slots: as the loci of the same particles in the list
    const loci = Particle.Loci(slots, OrderedSet.ofSortedArray([2, 3, 41, 85]));
    const same = Particle.Loci(list, OrderedSet.ofSortedArray([2, 3, 11, 25]));
    const s1 = Particle.getBoundingSphere(loci);
    const s2 = Particle.getBoundingSphere(same);
    expect(s1.radius).toBeCloseTo(s2.radius, 4);
    expect(Vec3.distance(s1.center, s2.center)).toBeLessThan(1e-4);
    // a gap's slot adds nothing
    const withGap = Particle.Loci(slots, OrderedSet.ofSortedArray([2, 3, 41, 60, 85]));
    expect(Particle.getBoundingSphere(withGap).radius).toBeCloseTo(s1.radius, 4);
    // all its live slots: its boundary; a page held whole by its list's, its particles inside, a page in part alike
    const all = Particle.Loci(slots, OrderedSet.ofSortedArray([...range(0, 10), ...range(40, 50), ...range(80, 90)]));
    const sAll = Particle.getBoundingSphere(all);
    expect(sAll.radius).toBeCloseTo(b.radius, 4);
    // into a sphere given: it set, the cached boundary not
    const into = Sphere3D.create(Vec3.create(1e6, 0, 0), 1);
    expect(Particle.getBoundingSphere(all, into)).toBe(into);
    expect(into.radius).toBeCloseTo(b.radius, 4);
    expect(Particle.getBoundary(slots).sphere.radius).toBeCloseTo(b.radius, 4);
    expect(Particle.getBoundary(slots).sphere.center[0]).not.toBe(1e6);
    const some = OrderedSet.ofSortedArray([...range(40, 50), 81, 84]);
    const sSome = Particle.getBoundingSphere(Particle.Loci(slots, some));
    const p = Vec3();
    OrderedSet.forEach(some, (v) => {
      Particle.position(p, Particle.Location(slots, v));
      expect(Vec3.distance(p, sSome.center)).toBeLessThanOrEqual(sSome.radius + 1e-3);
    });
    // no larger than the whole's of the same pages
    expect(sSome.radius).toBeLessThanOrEqual(sAll.radius + 1e-3);
  });

  it("keeps a loci's slots across versions, those of a page let go dropped", () => {
    const v1 = ParticleSlots.create(list, 100, pages());
    const v2 = ParticleSlots.update(v1, [], [40]);
    const loci = Particle.Loci(v1, OrderedSet.ofSortedArray([2, 41, 85]));
    const remapped = Particle.remapLoci(loci, v2);
    expect(remapped.particles).toBe(v2);
    expect(OrderedSet.toArray(remapped.indices)).toEqual([2, 85]);
    expect(Particle.getLabel(Particle.Loci(v2, OrderedSet.ofSingleton(85)))).toBe('particle 25');
  });
});
