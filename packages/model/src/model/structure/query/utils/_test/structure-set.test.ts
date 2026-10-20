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
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { SymmetryOperator } from '@molstar/core/math/geometry/symmetry-operator';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, StructureProperties as P, type StructureElement, Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import { Queries as Q, StructureQuery, StructureSelection } from '@molstar/model/model/structure/query';
import type { ElementIndex } from '../../../model/indexing.js';
import { StructureSymmetry } from '../../../structure/symmetry.js';
import { glycanCif } from '../../../structure/carbohydrates/_test/glycan.js';
import { structureUnion } from '../structure-set.js';

// Pins what `structureUnion` gives: the units (ids, invariant and chain-group ids, operators, element sets, whether a
// unit is the source's own), the symmetry groups and the root, for singletons, residues, overlapping and disjoint
// inputs, multi-unit and symmetry structures. The values were recorded before the 05-model-structure P1 change
// (.v6/perf-audit/05-model-structure.md), which must keep them; `foldUnion` is the implementation it replaced.

const examples = path.resolve(import.meta.dirname, '../../../../../../../../data/examples');

async function loadCif(name: string): Promise<Model> {
  const parsed = await CIF.parse(fs.readFileSync(path.join(examples, name), 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
}

/** The union before P1: one `SortedArray.union` per input unit. */
function foldUnion(source: Structure, structures: Structure[]) {
  if (structures.length === 0) return Structure.Empty;
  if (structures.length === 1) return structures[0];
  const unitMap = new Map<number, StructureElement.Set>();
  const fullUnits = new Set<number>();
  for (const { units } of structures) {
    for (const u of units) {
      if (unitMap.has(u.id)) {
        if (fullUnits.has(u.id)) continue;
        const merged = SortedArray.union(unitMap.get(u.id)!, u.elements);
        unitMap.set(u.id, merged);
        if (merged.length === source.unitMap.get(u.id).elements.length) fullUnits.add(u.id);
      } else {
        unitMap.set(u.id, u.elements);
        if (u.elements.length === source.unitMap.get(u.id).elements.length) fullUnits.add(u.id);
      }
    }
  }
  const builder = source.subsetBuilder(true);
  unitMap.forEach((elements, id) => builder.setUnit(id, elements));
  return builder.getStructure();
}

function hash(xs: ArrayLike<number>) {
  return (hashFnv32a(xs) >>> 0).toString(16);
}

/** Units and symmetry groups of a union of parts of `source`. */
function signature(source: Structure, s: Structure) {
  const rows = s.units.map((u) => [
    u.id,
    u.invariantId,
    u.chainGroupId,
    u.elements.length,
    hashFnv32a(u.elements),
    u === source.unitMap.get(u.id) ? 1 : 0,
  ]);
  const ops = s.units.map((u) => u.conformation.operator.name).join(',');
  const groups = s.unitSymmetryGroups.map((g) => g.units.map((u) => u.id));
  return [
    `units=${s.units.length} elements=${s.elementCount} sourceUnits=${rows.filter((r) => r[5]).length}`,
    `rows=${hash(rows.flat())} ops=${hash(Array.from(ops, (c) => c.charCodeAt(0)))}`,
    `groups=${groups.length}:${hash(groups.flatMap((g) => [g.length, ...g]))} root=${s.isEmpty || s.root === source.root}`,
  ].join(' ');
}

/** `a` and `b` are the same union: unit ids in order, elements, unit props, source units and symmetry groups. */
function expectSame(source: Structure, a: Structure, b: Structure) {
  expect(a.units.map((u) => u.id)).toEqual(b.units.map((u) => u.id));
  for (let i = 0; i < a.units.length; i++) {
    const u = a.units[i],
      v = b.units[i];
    expect(SortedArray.areEqual(u.elements, v.elements)).toBe(true);
    expect([u.invariantId, u.chainGroupId, u.kind, u.conformation.operator.name]).toEqual([
      v.invariantId,
      v.chainGroupId,
      v.kind,
      v.conformation.operator.name,
    ]);
    expect(u === source.unitMap.get(u.id)).toBe(v === source.unitMap.get(v.id));
  }
  const groups = (s: Structure) => s.unitSymmetryGroups.map((g) => g.units.map((u) => u.id));
  expect(groups(a)).toEqual(groups(b));
  expect(a.root === b.root).toBe(true);
}

/** A structure of `source` with the `pick(u)` indices of each unit that it has (skipping units where none is left). */
function part(source: Structure, pick: (u: Unit, k: number) => number[]) {
  const b = source.subsetBuilder(true);
  source.units.forEach((u, k) => {
    const idx = pick(u, k).filter((i) => i < u.elements.length);
    if (!idx.length) return;
    b.setUnit(u.id, SortedArray.ofSortedArray<ElementIndex>(Int32Array.from(idx, (i) => u.elements[i])));
  });
  return b.getStructure();
}

function range(start: number, end: number, step = 1) {
  const ret: number[] = [];
  for (let i = start; i < end; i += step) ret.push(i);
  return ret;
}

/** One structure per element: every `step`-th element of every unit, in a scrambled order, some twice. */
function singletons(source: Structure, step: number) {
  const ret: Structure[] = [];
  for (const u of source.units) {
    for (let i = 0; i < u.elements.length; i += step) {
      ret.push(Structure.create([u.getChild(SortedArray.ofSingleton(u.elements[i]))], { parent: source }));
    }
  }
  for (let i = ret.length - 1; i > 0; i--) {
    const j = (i * 7919) % (i + 1);
    [ret[i], ret[j]] = [ret[j], ret[i]];
  }
  for (let i = 0; i < ret.length; i += 5) ret.push(ret[i]);
  return ret;
}

/** All of a model's atoms in one unit, as MD and PDB files without chain ids give. */
function oneUnit(m: Model) {
  const all = SortedArray.ofSortedArray<ElementIndex>(
    Int32Array.from({ length: m.atomicHierarchy.atoms._rowCount }, (_, i) => i),
  );
  return Structure.create([
    Unit.create(0, 0, 0, Unit.Trait.MultiChain, Unit.Kind.Atomic, m, SymmetryOperator.Default, all),
  ]);
}

/** One structure per residue (every `step`-th by `label_seq_id`, waters and ligands included for 1). */
function residues(source: Structure, step = 1) {
  const residueTest: Parameters<typeof Q.generators.atoms>[0]['residueTest'] = (ctx) =>
    P.residue.label_seq_id(ctx.element) % step === 0;
  return StructureQuery.run(
    Q.generators.atoms({ residueTest, groupBy: (ctx) => P.residue.key(ctx.element) }),
    source,
  ) as StructureSelection.Sequence;
}

/** A seeded linear congruential generator in [0, 1). */
function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

describe('structureUnion', async () => {
  const m7qpd = await loadCif('7qpd.fw2.cif');
  const m1cbs = await loadCif('1cbs_updated.cif');
  const s7qpd = Structure.ofModel(m7qpd);
  const cell1cbs = await StructureSymmetry.buildSymmetryRange(
    Structure.ofModel(m1cbs),
    Vec3.create(0, 0, 0),
    Vec3.create(1, 0, 0),
  ).run();
  const one7qpd = oneUnit(m7qpd);

  type Case = [name: string, source: Structure, inputs: () => Structure[]];
  const cases: Case[] = [
    ['7qpd singletons', s7qpd, () => singletons(s7qpd, 3)],
    ['7qpd every singleton', s7qpd, () => singletons(s7qpd, 1)],
    ['7qpd one unit singletons', one7qpd, () => singletons(one7qpd, 2)],
    ['1cbs cells singletons', cell1cbs, () => singletons(cell1cbs, 4)],
    ['7qpd residues', s7qpd, () => residues(s7qpd).structures],
    ['7qpd every third residue', s7qpd, () => residues(s7qpd, 3).structures],
    ['7qpd one unit residues', one7qpd, () => residues(one7qpd).structures],
    ['7qpd one unit every third residue', one7qpd, () => residues(one7qpd, 3).structures],
    ['1cbs cells every second residue', cell1cbs, () => residues(cell1cbs, 2).structures],
    [
      '7qpd overlapping',
      s7qpd,
      () => [
        part(s7qpd, (u) => range(0, Math.ceil(u.elements.length * 0.6))),
        part(s7qpd, (u) => range(Math.floor(u.elements.length * 0.4), u.elements.length)),
        part(s7qpd, (u) => range(0, u.elements.length, 5)),
      ],
    ],
    [
      '7qpd overlapping partial',
      s7qpd,
      () => [
        part(s7qpd, (u) => range(0, u.elements.length, 3)),
        part(s7qpd, (u) => range(0, u.elements.length, 5)),
        part(s7qpd, (u) => range(10, Math.floor(u.elements.length / 2), 2)),
      ],
    ],
    [
      '7qpd disjoint, out of order',
      s7qpd,
      () => [
        part(s7qpd, (u, k) => (k % 2 ? range(Math.floor(u.elements.length / 2), u.elements.length) : [])),
        part(s7qpd, (u, k) => (k % 2 ? range(0, Math.floor(u.elements.length / 4)) : range(0, 7))),
        part(s7qpd, (u, k) => (k % 2 ? [] : range(9, 30, 3))),
      ],
    ],
    [
      '7qpd whole units and parts',
      s7qpd,
      () => [
        part(s7qpd, (u) => range(0, 10)),
        part(s7qpd, (u, k) => (k % 3 === 0 ? range(0, u.elements.length) : [])),
        part(s7qpd, (u) => range(5, 20)),
        part(s7qpd, (u, k) => (k % 3 === 1 ? range(20, u.elements.length) : [])),
      ],
    ],
    [
      '1cbs cells, same parts of every copy',
      cell1cbs,
      () => [
        part(cell1cbs, (u) => range(0, u.elements.length, 2)),
        part(cell1cbs, (u) => range(0, 40)),
        part(cell1cbs, (u) => range(1, u.elements.length, 6)),
      ],
    ],
    [
      '1cbs cells, different parts of each copy',
      cell1cbs,
      () => [
        part(cell1cbs, (u, k) => range(k % 4, u.elements.length, 3)),
        part(cell1cbs, (u, k) => (k % 2 ? range(0, u.elements.length) : [])),
        part(cell1cbs, (u) => range(0, 12)),
      ],
    ],
    [
      '7qpd touching windows',
      s7qpd,
      () => [part(s7qpd, () => range(0, 11)), part(s7qpd, () => range(10, 21)), part(s7qpd, () => range(20, 31))],
    ],
    [
      '7qpd touching windows, out of order',
      s7qpd,
      () => [part(s7qpd, () => range(20, 31)), part(s7qpd, () => range(0, 11)), part(s7qpd, () => range(10, 21))],
    ],
    [
      '7qpd one unit, few elements out of order',
      one7qpd,
      () => [
        part(one7qpd, () => range(500, 520)),
        part(one7qpd, () => range(100, 110)),
        part(one7qpd, () => range(105, 130)),
      ],
    ],
    [
      '7qpd one unit, large sets overlapping',
      one7qpd,
      () =>
        range(0, 10).map((k) => part(one7qpd, (u) => range(k * 100, k * 100 + Math.floor(u.elements.length * 0.9)))),
    ],
    ['7qpd two structures', s7qpd, () => [part(s7qpd, (u) => range(0, 50)), part(s7qpd, (u) => range(25, 90))]],
  ];

  it.each(cases)('matches the fold of unions: %s', (_, source, inputs) => {
    const xs = inputs();
    expectSame(source, structureUnion(source, xs), foldUnion(source, xs));
  });

  it('pins the unions', () => {
    const sig = Object.fromEntries(
      cases.map(([name, source, inputs]) => [name, signature(source, structureUnion(source, inputs()))]),
    );
    expect(sig).toMatchInlineSnapshot(`
      {
        "1cbs cells every second residue": "units=24 elements=5216 sourceUnits=16 rows=2db3cf35 ops=d6c8b9b9 groups=3:46ce97f7 root=true",
        "1cbs cells singletons": "units=24 elements=2432 sourceUnits=0 rows=43844d5d ops=d6c8b9b9 groups=3:46ce97f7 root=true",
        "1cbs cells, different parts of each copy": "units=24 elements=6566 sourceUnits=12 rows=6fb66079 ops=d6c8b9b9 groups=9:35cf313b root=true",
        "1cbs cells, same parts of every copy": "units=24 elements=6744 sourceUnits=8 rows=401c15d ops=d6c8b9b9 groups=3:46ce97f7 root=true",
        "7qpd disjoint, out of order": "units=7 elements=3594 sourceUnits=0 rows=2116633 ops=ffd43aee groups=7:7c4d33d9 root=true",
        "7qpd every singleton": "units=7 elements=11521 sourceUnits=7 rows=57d468e9 ops=ffd43aee groups=7:7c4d33d9 root=true",
        "7qpd every third residue": "units=7 elements=3945 sourceUnits=2 rows=42b99a0c ops=ffd43aee groups=7:7c4d33d9 root=true",
        "7qpd one unit every third residue": "units=1 elements=3945 sourceUnits=0 rows=cdafd2cb ops=93c84ec2 groups=1:eb741d64 root=true",
        "7qpd one unit residues": "units=1 elements=11521 sourceUnits=1 rows=2f1af5d2 ops=93c84ec2 groups=1:eb741d64 root=true",
        "7qpd one unit singletons": "units=1 elements=5761 sourceUnits=0 rows=d2503ce5 ops=93c84ec2 groups=1:eb741d64 root=true",
        "7qpd one unit, few elements out of order": "units=1 elements=50 sourceUnits=0 rows=c7ae1d41 ops=93c84ec2 groups=1:eb741d64 root=true",
        "7qpd one unit, large sets overlapping": "units=1 elements=11268 sourceUnits=0 rows=99432038 ops=93c84ec2 groups=1:eb741d64 root=true",
        "7qpd overlapping": "units=7 elements=11521 sourceUnits=7 rows=57d468e9 ops=ffd43aee groups=7:7c4d33d9 root=true",
        "7qpd overlapping partial": "units=7 elements=6891 sourceUnits=0 rows=1e304ef7 ops=ffd43aee groups=7:7c4d33d9 root=true",
        "7qpd residues": "units=7 elements=11521 sourceUnits=7 rows=57d468e9 ops=ffd43aee groups=7:7c4d33d9 root=true",
        "7qpd singletons": "units=7 elements=3843 sourceUnits=0 rows=ed508bf5 ops=ffd43aee groups=7:7c4d33d9 root=true",
        "7qpd touching windows": "units=7 elements=217 sourceUnits=0 rows=4b3346e5 ops=ffd43aee groups=7:7c4d33d9 root=true",
        "7qpd touching windows, out of order": "units=7 elements=217 sourceUnits=0 rows=4b3346e5 ops=ffd43aee groups=7:7c4d33d9 root=true",
        "7qpd two structures": "units=7 elements=572 sourceUnits=2 rows=3d0bbdff ops=ffd43aee groups=7:7c4d33d9 root=true",
        "7qpd whole units and parts": "units=7 elements=8157 sourceUnits=5 rows=61cdb6d9 ops=ffd43aee groups=7:7c4d33d9 root=true",
      }
    `);
  });

  it('returns an only input and the empty structure as they are', () => {
    const a = part(s7qpd, (u) => range(0, 3));
    expect(structureUnion(s7qpd, [a])).toBe(a);
    expect(structureUnion(s7qpd, [])).toBe(Structure.Empty);
  });

  it('matches the fold of unions on random inputs', () => {
    const rnd = random(42);
    for (const source of [s7qpd, cell1cbs, one7qpd]) {
      for (let trial = 0; trial < 40; trial++) {
        const n = 2 + Math.floor(rnd() * 6);
        const xs: Structure[] = [];
        for (let k = 0; k < n; k++) {
          const density = rnd();
          const kind = rnd();
          const s = part(source, (u) => {
            if (rnd() < 0.4) return [];
            const len = u.elements.length;
            if (kind < 0.15) return range(0, len);
            if (kind < 0.4) {
              const start = Math.floor(rnd() * len);
              return range(start, Math.min(len, start + 1 + Math.floor(rnd() * 40)));
            }
            return range(0, len).filter(() => rnd() < density * 0.3);
          });
          if (!s.isEmpty) xs.push(s);
        }
        expectSame(source, structureUnion(source, xs), foldUnion(source, xs));
      }
    }
  });

  it('keeps the callers: selections, queries and loci', () => {
    const q = (query: StructureQuery, s: Structure) => StructureQuery.run(query, s);
    const sig = (s: Structure, sel: StructureSelection) => signature(s, StructureSelection.unionStructure(sel));
    const residueSel = residues(s7qpd);
    const residue50 = Q.generators.atoms({ residueTest: (ctx) => P.residue.label_seq_id(ctx.element) === 50 });
    const within = q(Q.filters.within({ query: Q.generators.all, target: residue50, maxRadius: 5 }), s7qpd);
    const surroundings = q(Q.modifiers.includeSurroundings(residue50, { radius: 5, wholeResidues: true }), s7qpd);
    const loci = StructureSelection.toLociWithSourceUnits(residues(s7qpd, 3));
    expect({
      residues: sig(s7qpd, residueSel),
      union: sig(
        s7qpd,
        q(Q.modifiers.union(Q.generators.atoms({ groupBy: (ctx) => P.residue.key(ctx.element) })), s7qpd),
      ),
      within: `${StructureSelection.isSingleton(within)} ${sig(s7qpd, within)}`,
      surroundings: sig(s7qpd, surroundings),
      loci: loci.elements.map((e) => `${e.unit.id}:${OrderedSet.size(e.indices)}`).join(','),
    }).toMatchInlineSnapshot(`
      {
        "loci": "0:279,1:865,2:1107,3:734,4:838,5:39,6:83",
        "residues": "units=7 elements=11521 sourceUnits=7 rows=57d468e9 ops=ffd43aee groups=7:7c4d33d9 root=true",
        "surroundings": "units=5 elements=395 sourceUnits=0 rows=a287f0a3 ops=e153cf4a groups=5:e4f9d63c root=true",
        "union": "units=7 elements=11521 sourceUnits=7 rows=57d468e9 ops=ffd43aee groups=7:7c4d33d9 root=true",
        "within": "true units=5 elements=416 sourceUnits=0 rows=db776982 ops=e153cf4a groups=5:e4f9d63c root=true",
      }
    `);
  });

  // scaling, not a pin: the per-input merge copied the accumulated set for every singleton (quadratic, 16× for 4×)
  it('singletons in one unit scale near-linearly', { tags: ['heavy'] }, async () => {
    // one chain of 48k atoms
    const parsed = await CIF.parse(glycanCif({ strands: 40, length: 100, branchEvery: 3 })).run();
    if (parsed.isError) throw new Error(parsed.message);
    const source = Structure.ofModel(
      (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative,
    );
    const u = source.units[0];
    const take = (n: number) => {
      const step = Math.floor(u.elements.length / n);
      const xs: Structure[] = [];
      for (let i = 0; i < n; i++) {
        const e = u.elements[((i * 7919) % n) * step];
        xs.push(Structure.create([u.getChild(SortedArray.ofSingleton(e))], { parent: source }));
      }
      return xs;
    };
    const time = (xs: Structure[]) => {
      let best = Infinity;
      for (let i = 0; i < 9; i++) {
        const t0 = performance.now();
        structureUnion(source, xs);
        best = Math.min(best, performance.now() - t0);
      }
      return best;
    };
    const small = take(4000),
      large = take(16000);
    expect(time(large) / time(small)).toBeLessThan(10);
  });
});
