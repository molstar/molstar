/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { OrderedSet, SortedArray } from '@molstar/core/data/int';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Sphere3D } from '@molstar/core/math/geometry';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { InterUnitGraph } from '@molstar/core/math/graph/inter-unit-graph';
import { CIF } from '@molstar/io/reader/cif';
import { parseMol2 } from '@molstar/io/reader/mol2/parser';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { trajectoryFromMol2 } from '@molstar/model/formats/structure/mol2';
import { Structure, StructureElement, Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import type { ElementIndex } from '../../../model/indexing.js';
import { StructureSymmetry } from '../symmetry.js';
import { InterUnitBonds, type InterUnitEdgeProps } from '../unit/bonds/data.js';
import { computeInterUnitBonds } from '../unit/bonds/inter-compute.js';

// Pins what building structures and child structures gives today: bonds (pairs, orders, flags, keys), unit groupings,
// element sets and boundaries. The values were recorded before the 05-model-structure P4/P3/P5/P2 changes
// (.v6/perf-audit/05-model-structure.md), which must keep them.

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');

async function loadCif(name: string): Promise<Model> {
  const parsed = await CIF.parse(fs.readFileSync(path.join(examples, name), 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
}

async function loadMol2(name: string): Promise<Model> {
  const parsed = await parseMol2(fs.readFileSync(path.join(examples, name), 'utf8'), name).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMol2(parsed.result).run()).representative;
}

/** Splits every unit (at residue boundaries if `byResidue`) into chunks of about `size` elements, so bonds cross units. */
function splitUnits(s: Structure, size: number, byResidue = true): Structure {
  const units: Unit[] = [];
  let id = 0;
  for (const u of s.units) {
    const { elements } = u;
    const residueIndex = byResidue && Unit.isAtomic(u) ? u.residueIndex : undefined;
    let chunk: number[] = [];
    const flush = () => {
      if (chunk.length === 0) return;
      const set = SortedArray.ofSortedArray<ElementIndex>(Int32Array.from(chunk) as any);
      units.push(Unit.create(id, id, id, u.traits, u.kind, u.model, u.conformation.operator, set));
      id++;
      chunk = [];
    };
    for (let i = 0; i < elements.length; i++) {
      const e = elements[i];
      if (chunk.length >= size && (!residueIndex || residueIndex[e] !== residueIndex[elements[i - 1]])) flush();
      chunk.push(e);
    }
    flush();
  }
  return Structure.create(units);
}

function cell(s: Structure) {
  return StructureSymmetry.buildSymmetryRange(s, Vec3.create(0, 0, 0), Vec3.create(1, 0, 0)).run();
}

function hash(xs: ArrayLike<number>) {
  return (hashFnv32a(xs) >>> 0).toString(16);
}

function compareRows(a: number[], b: number[]) {
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return a[i] - b[i];
  return 0;
}

/** Inter-unit bonds: count, hash in edge order and hash of the sorted rows. */
function interSignature(bonds: InterUnitBonds) {
  const rows = bonds.edges.map((e) => [e.unitA, e.indexA, e.unitB, e.indexB, e.props.order, e.props.flag, e.props.key]);
  const flagCounts = new Map<number, number>();
  for (const r of rows) flagCounts.set(r[5], (flagCounts.get(r[5]) ?? 0) + 1);
  const flags = [...flagCounts].sort((a, b) => a[0] - b[0]).map(([f, c]) => `${f}:${c}`);
  return { count: rows.length, flags, ordered: hash(rows.flat()), sorted: hash([...rows].sort(compareRows).flat()) };
}

/** Intra-unit bonds of every atomic unit: the sorted (unit, a < b, order, flags, key) rows. */
function intraSignature(s: Structure) {
  const rows: number[][] = [];
  for (const u of s.units) {
    if (!Unit.isAtomic(u)) continue;
    const { a, b, edgeProps } = u.bonds;
    for (let i = 0; i < a.length; i++) {
      if (a[i] > b[i]) continue;
      rows.push([u.id, a[i], b[i], edgeProps.order[i], edgeProps.flags[i], edgeProps.key[i]]);
    }
  }
  return { count: rows.length, sorted: hash(rows.sort(compareRows).flat()) };
}

/** Units (ids, invariant and chain-group ids, operators, element sets) and symmetry groups. */
function unitSignature(s: Structure) {
  const rows = s.units.map((u) => [
    u.id,
    u.invariantId,
    u.chainGroupId,
    u.kind,
    u.elements.length,
    hashFnv32a(u.elements),
  ]);
  const ops = s.units.map((u) => u.conformation.operator.name).join(',');
  const groups = s.unitSymmetryGroups.map((g) => g.units.map((u) => u.id));
  return {
    units: s.units.length,
    elements: s.elementCount,
    rows: hash(rows.flat()),
    ops: hash(Array.from(ops, (c) => c.charCodeAt(0))),
    groups: hash(groups.flatMap((g) => [g.length, ...g])),
    groupCount: groups.length,
  };
}

/** Every second element of every unit (or the first half of each unit if `contiguous`), as a loci of `s`. */
function halfLoci(s: Structure, contiguous = false) {
  return StructureElement.Loci(
    s,
    s.units.map((u) => {
      const n = Math.ceil(u.elements.length / 2);
      const indices = contiguous
        ? Int32Array.from({ length: n }, (_, i) => i)
        : Int32Array.from({ length: n }, (_, i) => 2 * i);
      return { unit: u, indices: OrderedSet.ofSortedArray(SortedArray.ofSortedArray(indices)) as any };
    }),
  );
}

/** A child from the subset builder: whole units for `whole`, every second element of the rest. */
function subsetChild(s: Structure, whole: (u: Unit) => boolean) {
  const b = s.subsetBuilder(true);
  for (const u of s.units) {
    if (whole(u)) {
      b.setUnit(u.id, u.elements);
      continue;
    }
    b.beginUnit(u.id);
    for (let i = 0; i < u.elements.length; i += 2) b.addElement(u.elements[i]);
    b.commitUnit();
  }
  return b.getStructure();
}

/** The first residue (or element) of every unit, from the subset builder. */
function firstResidues(s: Structure) {
  const b = s.subsetBuilder(true);
  for (const u of s.units) {
    const { elements } = u;
    b.beginUnit(u.id);
    for (let i = 0; i < elements.length; i++) {
      if (Unit.isAtomic(u) && i > 0 && u.residueIndex[elements[i]] !== u.residueIndex[elements[0]]) break;
      b.addElement(elements[i]);
      if (!Unit.isAtomic(u)) break;
    }
    b.commitUnit();
  }
  return b.getStructure();
}

/** The residues at the inter-unit bonds of `s` (computed aside, so `s` keeps no bonds), from the subset builder. */
function bondedResidues(s: Structure) {
  const byUnit = new Map<number, Set<number>>();
  for (const e of computeInterUnitBonds(s).edges) {
    for (const [id, index] of [
      [e.unitA, e.indexA],
      [e.unitB, e.indexB],
    ]) {
      const u = s.unitMap.get(id) as Unit.Atomic;
      if (!byUnit.has(id)) byUnit.set(id, new Set());
      byUnit.get(id)!.add(u.residueIndex[u.elements[index]]);
    }
  }
  const b = s.subsetBuilder(true);
  for (const u of s.units) {
    const residues = byUnit.get(u.id);
    if (!residues || !Unit.isAtomic(u)) continue;
    b.beginUnit(u.id);
    for (let i = 0; i < u.elements.length; i++)
      if (residues.has(u.residueIndex[u.elements[i]])) b.addElement(u.elements[i]);
    b.commitUnit();
  }
  return b.getStructure();
}

/** The largest distance by which an element (centre plus radius) of `s` lies outside `sphere`. */
function worstOutside(sphere: { center: Vec3; radius: number }, s: Structure) {
  let worst = -Infinity;
  const p = Vec3();
  for (const u of s.units) {
    for (let i = 0; i < u.elements.length; i++) {
      const e = u.elements[i];
      u.conformation.position(e, p);
      worst = Math.max(worst, Vec3.distance(p, sphere.center) + u.conformation.r(e) - sphere.radius);
    }
  }
  return worst;
}

function round(x: number, digits = 3) {
  const f = 10 ** digits;
  return Math.round(x * f) / f;
}

function sphereOf(s: { center: Vec3; radius: number }) {
  return [...Array.from(s.center, (x) => round(x)), round(s.radius)];
}

/** Recorded values: `PIN=<file> pnpm test:unit <this file>` writes the current ones to <file>. */
const PIN = !!process.env.PIN;
const pinned: Record<string, unknown> = {};
function check(key: string, actual: unknown, expected: Record<string, unknown>) {
  if (PIN) pinned[key] = actual;
  else expect(compact(actual), key).toEqual(expected[key]);
}

/** `{ a: 1, b: { c: 2 } }` → `'a=1 b=c=2'`, arrays comma-separated. */
function compact(v: unknown): string {
  if (Array.isArray(v)) return v.join(',');
  if (v && typeof v === 'object')
    return Object.entries(v)
      .map(([k, x]) => `${k}=${compact(x)}`)
      .join(' ');
  return String(v);
}

type Case = [name: string, root: () => Promise<Structure>];
type Child = [name: string, child: (root: Structure) => Structure];

describe('structure building', async () => {
  const m1bna = await loadCif('1bna_confal_pyramids.cif');
  const m7qpd = await loadCif('7qpd.fw2.cif');
  const m1cbs = await loadCif('1cbs_updated.cif');
  const m1grm = await loadCif('1grm_updated.cif');
  const mAce2 = await loadMol2('ace2-hit.mol2');

  const cases: Case[] = [
    // struct_conn hydrogen bonds between symmetry copies
    ['1bna cell', () => cell(Structure.ofModel(m1bna))],
    // struct_conn covalent and disulfide bonds between units
    ['7qpd', async () => Structure.ofModel(m7qpd)],
    // computed covalent bonds between units, under symmetry operators; ligand and waters
    ['1cbs split cell', () => cell(splitUnits(Structure.ofModel(m1cbs), 300))],
    // index-pair bonds between units
    ['ace2 split', async () => splitUnits(Structure.ofModel(mAce2), 12, false)],
    // struct_conn bonds inside units
    ['1grm', async () => Structure.ofModel(m1grm)],
  ];

  const children: Child[] = [
    ['half', (s) => subsetChild(s, () => false)],
    ['polymer whole, rest half', (s) => subsetChild(s, (u) => u.polymerElements.length > 0)],
    ['bundle half', (s) => StructureElement.Bundle.toStructure(StructureElement.Bundle.fromLoci(halfLoci(s)), s)],
    ['loci half', (s) => StructureElement.Loci.toStructure(halfLoci(s))],
    [
      'bundle first halves',
      (s) => StructureElement.Bundle.toStructure(StructureElement.Bundle.fromLoci(halfLoci(s, true)), s),
    ],
    ['loci first halves', (s) => StructureElement.Loci.toStructure(halfLoci(s, true))],
    // small children of large units: the subsets search the parent instead of walking it
    ['first residue of each unit', (s) => firstResidues(s)],
    ['residues at inter-unit bonds', (s) => bondedResidues(s)],
  ];

  it.each(cases)('root bonds: %s', async (name, root) => {
    const s = await root();
    check(`${name}: units`, unitSignature(s), expected);
    check(`${name}: inter`, interSignature(s.interUnitBonds), expected);
    check(`${name}: intra`, intraSignature(s), expected);
  });

  // child bonds are the same whether the root's bonds exist (subset) or not (computed)
  for (const rootBonds of [false, true]) {
    it.each(cases)(`child bonds, root bonds ${rootBonds ? 'computed' : 'absent'}: %s`, async (name, root) => {
      for (const [childName, child] of children) {
        const s = await root();
        if (rootBonds) {
          void s.interUnitBonds;
          for (const u of s.units) if (Unit.isAtomic(u)) void u.bonds;
        }
        const c = child(s);
        const key = `${name} / ${childName}`;
        check(`${key}: units`, unitSignature(c), expected);
        const { count, flags, sorted } = interSignature(c.interUnitBonds);
        check(`${key}: inter`, { count, flags, sorted }, expected);
        check(`${key}: intra`, intraSignature(c), expected);
      }
    });
  }

  describe('child inter-unit bonds', () => {
    /** The root's units with one made-up bond (that no computation would find) between the first atoms of units 0, 1. */
    async function rootWithMadeUpBond() {
      const s = Structure.ofModel(m7qpd);
      const builder = new InterUnitGraph.Builder<number, StructureElement.UnitIndex, InterUnitEdgeProps>();
      builder.startUnitPair(s.units[0].id, s.units[1].id);
      builder.add(0 as StructureElement.UnitIndex, 0 as StructureElement.UnitIndex, { order: 1, flag: 1, key: -1 });
      builder.finishUnitPair();
      return Structure.create(s.units, { interUnitBonds: new InterUnitBonds(builder.getMap()) });
    }

    it("restricts the parent's bonds", async () => {
      const root = await rootWithMadeUpBond();
      // every second element keeps index 0 of every unit
      const half = subsetChild(root, () => false);
      expect(half.interUnitBonds.edgeCount).toBe(2);
      expect(half.interUnitBonds.getEdge(0 as any, root.units[0].id, 0 as any, root.units[1].id)).toBeDefined();
      // whole units: the parent's bonds as they are
      expect(subsetChild(root, () => true).interUnitBonds).toBe(root.interUnitBonds);
      // without the first element of unit 0
      const b = root.subsetBuilder(true);
      for (const u of root.units) b.setUnit(u.id, u === root.units[0] ? u.elements.slice(1) : u.elements);
      expect(b.getStructure().interUnitBonds.edgeCount).toBe(0);
    });

    it('computes when the units are transformed', async () => {
      const root = await rootWithMadeUpBond();
      void root.interUnitBonds;
      // same unit ids, other operators: computed, so the made-up bond is gone and 7qpd's own 6 are there
      const moved = Structure.transform(root, Mat4.fromTranslation(Mat4(), Vec3.create(10, 0, 0)));
      expect(moved.parent).toBe(root);
      expect(interSignature(moved.interUnitBonds).sorted).toBe(
        interSignature(Structure.ofModel(m7qpd).interUnitBonds).sorted,
      );
    });
  });

  // symmetry copies of a child unit share its props (lookups, bonds, boundary) whichever way the child is built
  it.each(children)('child props shared across symmetry copies: %s', async (_name, child) => {
    for (const [, root] of [cases[0], cases[2]]) {
      const c = child(await root());
      expect(new Set(c.units.map((u) => u.props)).size).toBe(c.unitSymmetryGroups.length);
    }
  });

  it.each(cases)('boundaries: %s', async (name, root) => {
    const s = await root();
    // unit spheres: exact
    check(`${name}: unit spheres`, hash(s.units.flatMap((u) => sphereOf(u.boundary.sphere))), expected);
    // structure box and sphere: exact (these fixtures stay under the 10k points of the 98-direction helper)
    const { box, sphere } = s.boundary;
    check(`${name}: box`, [...Array.from(box.min, (x) => round(x)), ...Array.from(box.max, (x) => round(x))], expected);
    const [cx, cy, cz, r] = sphereOf(sphere);
    if (PIN) pinned[`${name}: sphere`] = [cx, cy, cz, r];
    else expect([cx, cy, cz, r]).toEqual(expected[`${name}: sphere`]);
    // loci sphere of every second element: exact, encloses every selected element
    const loci = halfLoci(s);
    const ls = StructureElement.Loci.getBoundingSphere(loci);
    const lr = round(ls.radius);
    if (PIN) pinned[`${name}: half loci radius`] = lr;
    else {
      expect(lr).toBe(expected[`${name}: half loci radius`]);
      expect(worstOutside(ls, StructureElement.Loci.toStructure(loci))).toBeLessThan(1e-3);
    }
  });

  // above 10k points the spheres use 14 directions: the box is the same, the radius within 3 % of the 98-direction one
  // and the atoms are enclosed as closely; the centre may move (11 Å, 7 % of the radius here), as both are
  // approximations of the smallest sphere (05-model-structure P2)
  it('boundaries of a large structure', async () => {
    const s = await StructureSymmetry.buildSymmetryRange(
      Structure.ofModel(m1cbs),
      Vec3.create(0, 0, 0),
      Vec3.create(2, 2, 2),
    ).run();
    const reference = new BoundaryHelper('98');
    const tmp = Sphere3D();
    const unitSphere = (u: Unit) => Sphere3D.transform(tmp, u.boundary.sphere, u.conformation.operator.matrix);
    reference.reset();
    for (const u of s.units) reference.includeSphere(unitSphere(u));
    reference.finishedIncludeStep();
    for (const u of s.units) reference.radiusSphere(unitSphere(u));
    const expectedSphere = reference.getSphere();

    const { box, sphere } = s.boundary;
    // the 14-direction path: at most 28 distinct extrema
    expect(new Set(sphere.extrema!.map((e) => e.join())).size).toBeLessThanOrEqual(28);
    const expectedBox = reference.getBox();
    for (let i = 0; i < 3; i++) {
      expect(box.min[i]).toBeCloseTo(expectedBox.min[i], 6);
      expect(box.max[i]).toBeCloseTo(expectedBox.max[i], 6);
    }
    expect(Math.abs(sphere.radius - expectedSphere.radius) / expectedSphere.radius).toBeLessThan(0.03);
    expect(worstOutside(sphere, s)).toBeLessThan(worstOutside(expectedSphere, s) + 0.01);

    const loci = halfLoci(s);
    expect(
      worstOutside(StructureElement.Loci.getBoundingSphere(loci), StructureElement.Loci.toStructure(loci)),
    ).toBeLessThan(1e-3);
  });

  if (PIN) {
    it('print', () => {
      fs.writeFileSync(process.env.PIN!, JSON.stringify(pinned, null, 2));
    });
  }
});

const expected: Record<string, unknown> = {
  '1bna cell: units': 'units=24 elements=4528 rows=926e64f5 ops=14036e05 groups=46ce97f7 groupCount=3',
  '1bna cell: inter': 'count=512 flags=4:512 ordered=bdc80e85 sorted=2ebef715',
  '1bna cell: intra': 'count=4352 sorted=37340355',
  '7qpd: units': 'units=7 elements=11521 rows=6dafd90e ops=ffd43aee groups=7c4d33d9 groupCount=7',
  '7qpd: inter': 'count=6 flags=1:4,9:2 ordered=9ec72f2b sorted=353437d1',
  '7qpd: intra': 'count=11854 sorted=6ab707d5',
  '1cbs split cell: units': 'units=48 elements=9704 rows=89c11605 ops=70a54e4f groups=cfa8213d groupCount=6',
  '1cbs split cell: inter': 'count=48 flags=33:48 ordered=e1b88ea5 sorted=5f95c8c5',
  '1cbs split cell: intra': 'count=9008 sorted=99b7c635',
  'ace2 split: units': 'units=4 elements=37 rows=578b6ba2 ops=cef4bbf3 groups=9f25f271 groupCount=4',
  'ace2 split: inter': 'count=20 flags=1:12,17:8 ordered=d76b2aef sorted=dffd54db',
  'ace2 split: intra': 'count=31 sorted=82125e7f',
  '1grm: units': 'units=2 elements=272 rows=635bf8b6 ops=61a2fd77 groups=206c2eda groupCount=2',
  '1grm: inter': 'count=0 flags= ordered=811c9dc5 sorted=811c9dc5',
  '1grm: intra': 'count=286 sorted=6623881c',
  '1bna cell / half: units': 'units=24 elements=2272 rows=b928d2c5 ops=14036e05 groups=46ce97f7 groupCount=3',
  '1bna cell / half: inter': 'count=128 flags=4:128 sorted=2f35025',
  '1bna cell / half: intra': 'count=704 sorted=d875d3c5',
  '1bna cell / polymer whole, rest half: units':
    'units=24 elements=4208 rows=32a25da5 ops=14036e05 groups=46ce97f7 groupCount=3',
  '1bna cell / polymer whole, rest half: inter': 'count=512 flags=4:512 sorted=2ebef715',
  '1bna cell / polymer whole, rest half: intra': 'count=4352 sorted=37340355',
  '1bna cell / bundle half: units': 'units=24 elements=2272 rows=b928d2c5 ops=14036e05 groups=46ce97f7 groupCount=3',
  '1bna cell / bundle half: inter': 'count=128 flags=4:128 sorted=2f35025',
  '1bna cell / bundle half: intra': 'count=704 sorted=d875d3c5',
  '1bna cell / loci half: units': 'units=24 elements=2272 rows=b928d2c5 ops=14036e05 groups=46ce97f7 groupCount=3',
  '1bna cell / loci half: inter': 'count=128 flags=4:128 sorted=2f35025',
  '1bna cell / loci half: intra': 'count=704 sorted=d875d3c5',
  '1bna cell / bundle first halves: units':
    'units=24 elements=2272 rows=ba3a9a8d ops=14036e05 groups=46ce97f7 groupCount=3',
  '1bna cell / bundle first halves: inter': 'count=0 flags= sorted=811c9dc5',
  '1bna cell / bundle first halves: intra': 'count=2192 sorted=b739f7fd',
  '1bna cell / loci first halves: units':
    'units=24 elements=2272 rows=ba3a9a8d ops=14036e05 groups=46ce97f7 groupCount=3',
  '1bna cell / loci first halves: inter': 'count=0 flags= sorted=811c9dc5',
  '1bna cell / loci first halves: intra': 'count=2192 sorted=b739f7fd',
  '7qpd / half: units': 'units=7 elements=5763 rows=d51cab90 ops=ffd43aee groups=7c4d33d9 groupCount=7',
  '7qpd / half: inter': 'count=2 flags=1:2 sorted=8c245ef7',
  '7qpd / half: intra': 'count=1094 sorted=f0c3a3d7',
  '7qpd / polymer whole, rest half: units':
    'units=7 elements=11461 rows=4c56f883 ops=ffd43aee groups=7c4d33d9 groupCount=7',
  '7qpd / polymer whole, rest half: inter': 'count=6 flags=1:4,9:2 sorted=353437d1',
  '7qpd / polymer whole, rest half: intra': 'count=11733 sorted=f388f444',
  '7qpd / bundle half: units': 'units=7 elements=5763 rows=d51cab90 ops=ffd43aee groups=7c4d33d9 groupCount=7',
  '7qpd / bundle half: inter': 'count=2 flags=1:2 sorted=8c245ef7',
  '7qpd / bundle half: intra': 'count=1094 sorted=f0c3a3d7',
  '7qpd / loci half: units': 'units=7 elements=5763 rows=d51cab90 ops=ffd43aee groups=7c4d33d9 groupCount=7',
  '7qpd / loci half: inter': 'count=2 flags=1:2 sorted=8c245ef7',
  '7qpd / loci half: intra': 'count=1094 sorted=f0c3a3d7',
  '7qpd / bundle first halves: units': 'units=7 elements=5763 rows=aac0c7b6 ops=ffd43aee groups=7c4d33d9 groupCount=7',
  '7qpd / bundle first halves: inter': 'count=4 flags=1:2,9:2 sorted=d461e6f',
  '7qpd / bundle first halves: intra': 'count=5917 sorted=8d4e8153',
  '7qpd / loci first halves: units': 'units=7 elements=5763 rows=aac0c7b6 ops=ffd43aee groups=7c4d33d9 groupCount=7',
  '7qpd / loci first halves: inter': 'count=4 flags=1:2,9:2 sorted=d461e6f',
  '7qpd / loci first halves: intra': 'count=5917 sorted=8d4e8153',
  '1cbs split cell / half: units': 'units=48 elements=4864 rows=681f49c5 ops=70a54e4f groups=cfa8213d groupCount=6',
  '1cbs split cell / half: inter': 'count=32 flags=33:32 sorted=d4b5cd45',
  '1cbs split cell / half: intra': 'count=848 sorted=b4f1c6b5',
  '1cbs split cell / polymer whole, rest half: units':
    'units=48 elements=9216 rows=61eb5d5 ops=70a54e4f groups=cfa8213d groupCount=6',
  '1cbs split cell / polymer whole, rest half: inter': 'count=48 flags=33:48 sorted=5f95c8c5',
  '1cbs split cell / polymer whole, rest half: intra': 'count=8856 sorted=4c460b25',
  '1cbs split cell / bundle half: units':
    'units=48 elements=4864 rows=681f49c5 ops=70a54e4f groups=cfa8213d groupCount=6',
  '1cbs split cell / bundle half: inter': 'count=32 flags=33:32 sorted=d4b5cd45',
  '1cbs split cell / bundle half: intra': 'count=848 sorted=b4f1c6b5',
  '1cbs split cell / loci half: units':
    'units=48 elements=4864 rows=681f49c5 ops=70a54e4f groups=cfa8213d groupCount=6',
  '1cbs split cell / loci half: inter': 'count=32 flags=33:32 sorted=d4b5cd45',
  '1cbs split cell / loci half: intra': 'count=848 sorted=b4f1c6b5',
  '1cbs split cell / bundle first halves: units':
    'units=48 elements=4864 rows=c091c6a5 ops=70a54e4f groups=cfa8213d groupCount=6',
  '1cbs split cell / bundle first halves: inter': 'count=0 flags= sorted=811c9dc5',
  '1cbs split cell / bundle first halves: intra': 'count=4512 sorted=69182535',
  '1cbs split cell / loci first halves: units':
    'units=48 elements=4864 rows=c091c6a5 ops=70a54e4f groups=cfa8213d groupCount=6',
  '1cbs split cell / loci first halves: inter': 'count=0 flags= sorted=811c9dc5',
  '1cbs split cell / loci first halves: intra': 'count=4512 sorted=69182535',
  'ace2 split / half: units': 'units=4 elements=19 rows=2dce34a2 ops=cef4bbf3 groups=9f25f271 groupCount=4',
  'ace2 split / half: inter': 'count=6 flags=1:4,17:2 sorted=e58f3c8f',
  'ace2 split / half: intra': 'count=0 sorted=811c9dc5',
  'ace2 split / polymer whole, rest half: units':
    'units=4 elements=19 rows=2dce34a2 ops=cef4bbf3 groups=9f25f271 groupCount=4',
  'ace2 split / polymer whole, rest half: inter': 'count=6 flags=1:4,17:2 sorted=e58f3c8f',
  'ace2 split / polymer whole, rest half: intra': 'count=0 sorted=811c9dc5',
  'ace2 split / bundle half: units': 'units=4 elements=19 rows=2dce34a2 ops=cef4bbf3 groups=9f25f271 groupCount=4',
  'ace2 split / bundle half: inter': 'count=6 flags=1:4,17:2 sorted=e58f3c8f',
  'ace2 split / bundle half: intra': 'count=0 sorted=811c9dc5',
  'ace2 split / loci half: units': 'units=4 elements=19 rows=2dce34a2 ops=cef4bbf3 groups=9f25f271 groupCount=4',
  'ace2 split / loci half: inter': 'count=6 flags=1:4,17:2 sorted=e58f3c8f',
  'ace2 split / loci half: intra': 'count=0 sorted=811c9dc5',
  'ace2 split / bundle first halves: units':
    'units=4 elements=19 rows=4c9fb4c5 ops=cef4bbf3 groups=9f25f271 groupCount=4',
  'ace2 split / bundle first halves: inter': 'count=8 flags=1:6,17:2 sorted=36e58423',
  'ace2 split / bundle first halves: intra': 'count=12 sorted=96cb31b',
  'ace2 split / loci first halves: units':
    'units=4 elements=19 rows=4c9fb4c5 ops=cef4bbf3 groups=9f25f271 groupCount=4',
  'ace2 split / loci first halves: inter': 'count=8 flags=1:6,17:2 sorted=36e58423',
  'ace2 split / loci first halves: intra': 'count=12 sorted=96cb31b',
  '1grm / half: units': 'units=2 elements=136 rows=57a23dde ops=61a2fd77 groups=206c2eda groupCount=2',
  '1grm / half: inter': 'count=0 flags= sorted=811c9dc5',
  '1grm / half: intra': 'count=46 sorted=1c647126',
  '1grm / polymer whole, rest half: units':
    'units=2 elements=272 rows=635bf8b6 ops=61a2fd77 groups=206c2eda groupCount=2',
  '1grm / polymer whole, rest half: inter': 'count=0 flags= sorted=811c9dc5',
  '1grm / polymer whole, rest half: intra': 'count=286 sorted=6623881c',
  '1grm / bundle half: units': 'units=2 elements=136 rows=57a23dde ops=61a2fd77 groups=206c2eda groupCount=2',
  '1grm / bundle half: inter': 'count=0 flags= sorted=811c9dc5',
  '1grm / bundle half: intra': 'count=46 sorted=1c647126',
  '1grm / loci half: units': 'units=2 elements=136 rows=57a23dde ops=61a2fd77 groups=206c2eda groupCount=2',
  '1grm / loci half: inter': 'count=0 flags= sorted=811c9dc5',
  '1grm / loci half: intra': 'count=46 sorted=1c647126',
  '1grm / bundle first halves: units': 'units=2 elements=136 rows=5a061b76 ops=61a2fd77 groups=206c2eda groupCount=2',
  '1grm / bundle first halves: inter': 'count=0 flags= sorted=811c9dc5',
  '1grm / bundle first halves: intra': 'count=138 sorted=df350ae0',
  '1grm / loci first halves: units': 'units=2 elements=136 rows=5a061b76 ops=61a2fd77 groups=206c2eda groupCount=2',
  '1grm / loci first halves: inter': 'count=0 flags= sorted=811c9dc5',
  '1grm / loci first halves: intra': 'count=138 sorted=df350ae0',
  '1bna cell: unit spheres': 'fc0869f5',
  '1bna cell: box': '-11.333,-15.52,-13.081,61.073,55.91,79.281',
  '1bna cell: sphere': [19.838, 18.015, 34.215, 54.83],
  '1bna cell: half loci radius': 53.965,
  '7qpd: unit spheres': '365ed797',
  '7qpd: box': '147.773,119.237,144.23,238.404,247.235,266.855',
  '7qpd: sphere': [192.432, 186.805, 205.346, 73.882],
  '7qpd: half loci radius': 73.144,
  '1cbs split cell: unit spheres': 'a73e2fd5',
  '1cbs split cell: box': '-15.796,-12.661,-12.216,107.096,60.221,89.826',
  '1cbs split cell: sphere': [43.928, 26.527, 40.729, 71.729],
  '1cbs split cell: half loci radius': 71.334,
  'ace2 split: unit spheres': 'ab5765b2',
  'ace2 split: box': '62.049,76.941,27.882,71.887,90.831,40.835',
  'ace2 split: sphere': [67.064, 84.401, 34.51, 7.749],
  'ace2 split: half loci radius': 7.895,
  '1grm: unit spheres': 'c1defe0f',
  '1grm: box': '-10.448,-14.779,-6.382,10.421,14.756,14.354',
  '1grm: sphere': [-0.022, 0.093, 3.901, 16.73],
  '1grm: half loci radius': 16.72,
  '1bna cell / first residue of each unit: units':
    'units=24 elements=264 rows=37017ead ops=14036e05 groups=46ce97f7 groupCount=3',
  '1bna cell / first residue of each unit: inter': 'count=0 flags= sorted=811c9dc5',
  '1bna cell / first residue of each unit: intra': 'count=272 sorted=7d83b0bd',
  '1bna cell / residues at inter-unit bonds: units':
    'units=16 elements=3888 rows=142a192d ops=cdaef8db groups=43ea81a5 groupCount=2',
  '1bna cell / residues at inter-unit bonds: inter': 'count=512 flags=4:512 sorted=2ebef715',
  '1bna cell / residues at inter-unit bonds: intra': 'count=4352 sorted=37340355',
  '7qpd / first residue of each unit: units':
    'units=7 elements=53 rows=1c48e514 ops=ffd43aee groups=7c4d33d9 groupCount=7',
  '7qpd / first residue of each unit: inter': 'count=0 flags= sorted=811c9dc5',
  '7qpd / first residue of each unit: intra': 'count=48 sorted=23f4e1b3',
  '7qpd / residues at inter-unit bonds: units':
    'units=5 elements=56 rows=6589f47f ops=e153cf4a groups=a4438b3d groupCount=5',
  '7qpd / residues at inter-unit bonds: inter': 'count=6 flags=1:4,9:2 sorted=580e2015',
  '7qpd / residues at inter-unit bonds: intra': 'count=52 sorted=d8625be4',
  '1cbs split cell / first residue of each unit: units':
    'units=48 elements=408 rows=9e3895e5 ops=70a54e4f groups=cfa8213d groupCount=6',
  '1cbs split cell / first residue of each unit: inter': 'count=0 flags= sorted=811c9dc5',
  '1cbs split cell / first residue of each unit: intra': 'count=376 sorted=3f1dc565',
  '1cbs split cell / residues at inter-unit bonds: units':
    'units=32 elements=336 rows=34364175 ops=fc9023d3 groups=fa776e5 groupCount=4',
  '1cbs split cell / residues at inter-unit bonds: inter': 'count=48 flags=33:48 sorted=39be0a5',
  '1cbs split cell / residues at inter-unit bonds: intra': 'count=296 sorted=78871925',
  'ace2 split / first residue of each unit: units':
    'units=4 elements=37 rows=578b6ba2 ops=cef4bbf3 groups=9f25f271 groupCount=4',
  'ace2 split / first residue of each unit: inter': 'count=20 flags=1:12,17:8 sorted=dffd54db',
  'ace2 split / first residue of each unit: intra': 'count=31 sorted=82125e7f',
  'ace2 split / residues at inter-unit bonds: units':
    'units=4 elements=37 rows=578b6ba2 ops=cef4bbf3 groups=9f25f271 groupCount=4',
  'ace2 split / residues at inter-unit bonds: inter': 'count=20 flags=1:12,17:8 sorted=dffd54db',
  'ace2 split / residues at inter-unit bonds: intra': 'count=31 sorted=82125e7f',
  '1grm / first residue of each unit: units':
    'units=2 elements=18 rows=42c69092 ops=61a2fd77 groups=206c2eda groupCount=2',
  '1grm / first residue of each unit: inter': 'count=0 flags= sorted=811c9dc5',
  '1grm / first residue of each unit: intra': 'count=16 sorted=e3e71edd',
  '1grm / residues at inter-unit bonds: units':
    'units=0 elements=0 rows=811c9dc5 ops=811c9dc5 groups=811c9dc5 groupCount=0',
  '1grm / residues at inter-unit bonds: inter': 'count=0 flags= sorted=811c9dc5',
  '1grm / residues at inter-unit bonds: intra': 'count=0 sorted=811c9dc5',
};
