/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import * as fs from 'node:fs';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, Unit } from '@molstar/model/model/structure';
import { MVSAnnotation } from '../../components/annotation-prop.js';
import { substructureFromRows } from '@molstar/model/model/structure/structure/element/schema-selection/row-structure';
import type { Model } from '@molstar/model/model/structure';
import type { SchemaItem } from '@molstar/model/model/structure/structure/element/schema';
import type { MVSAnnotationRow } from '../schemas.js';
import { IdentityOperator, mulberry32, randInt } from './annotation-rows-fixture.js';

/** `data/examples/ihm-coarse-8zzh.cif`: 4 coarse chains of 882 elements (A, C: spheres; B, D: gaussians), 2 entities, no atoms. */
async function loadCoarseModel() {
  const parsed = await CIF.parseText(fs.readFileSync('data/examples/ihm-coarse-8zzh.cif', 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return trajectory.representative;
}

/** An independent definition of "a coarse element satisfies a row" (the one of `StructureElement.Schema`): the chain and the entity of the element are
 * as given, a residue or the bounds of a range overlap the sequence the element spans (`seq_id_begin` to `seq_id_end`, both inclusive), all the criteria
 * which are given hold, and a criterion which only atoms have (`auth_*`, the comp ids, atom ids, ...) is not satisfied by anything. A range which ends
 * before it begins selects nothing. Rows apply to the one operator of the structure (the identity). */
function createCoarseOracle(model: Model) {
  const AtomOnly = [
    'auth_asym_id',
    'auth_seq_id',
    'pdbx_PDB_ins_code',
    'beg_auth_seq_id',
    'end_auth_seq_id',
    'beg_auth_ins_code',
    'end_auth_ins_code',
    'label_comp_id',
    'auth_comp_id',
    'residue_index',
    'label_atom_id',
    'auth_atom_id',
    'type_symbol',
    'atom_id',
    'atom_index',
  ] as const;
  const defined = (x: unknown) => x !== undefined && x !== null;
  function satisfies(row: SchemaItem, e: Model['coarseHierarchy']['spheres'], i: number): boolean {
    if (AtomOnly.some((k) => defined(row[k]))) return false;
    if (defined(row.instance_id) && row.instance_id !== IdentityOperator.instanceId) return false;
    if (defined(row.operator_name) && row.operator_name !== IdentityOperator.name) return false;
    if (defined(row.label_asym_id) && row.label_asym_id !== e.asym_id.value(i)) return false;
    if (defined(row.label_entity_id) && row.label_entity_id !== e.entity_id.value(i)) return false;
    const begin = e.seq_id_begin.value(i);
    const end = e.seq_id_end.value(i);
    if (defined(row.label_seq_id) && !(begin <= row.label_seq_id && row.label_seq_id <= end)) return false;
    if (defined(row.beg_label_seq_id) && defined(row.end_label_seq_id) && row.beg_label_seq_id > row.end_label_seq_id)
      return false;
    if (defined(row.beg_label_seq_id) && !(end >= row.beg_label_seq_id)) return false;
    if (defined(row.end_label_seq_id) && !(begin <= row.end_label_seq_id)) return false;
    return true;
  }
  return {
    /** The index of the last of `rows` which each element satisfies, -1 for none */
    rowOfEach(rows: readonly SchemaItem[]) {
      const of = (e: Model['coarseHierarchy']['spheres']) => {
        const result = new Array<number>(e.count).fill(-1);
        rows.forEach((row, iRow) => {
          for (let i = 0; i < e.count; i++) if (satisfies(row, e, i)) result[i] = iRow;
        });
        return result;
      };
      return { spheres: of(model.coarseHierarchy.spheres), gaussians: of(model.coarseHierarchy.gaussians) };
    },
    /** The units (by kind: the elements) which satisfy at least one of `rows` */
    elementsOf(rows: readonly SchemaItem[]) {
      const rowOf = this.rowOfEach(rows);
      return {
        spheres: rowOf.spheres.flatMap((v, i) => (v >= 0 ? [i] : [])),
        gaussians: rowOf.gaussians.flatMap((v, i) => (v >= 0 ? [i] : [])),
      };
    },
  };
}
/** What `substructureFromRows` selected, as the elements of spheres and gaussians */
function selectedElements(structure: Structure) {
  const result = { spheres: [] as number[], gaussians: [] as number[] };
  for (const unit of structure.units) {
    const target = unit.kind === Unit.Kind.Spheres ? result.spheres : result.gaussians;
    for (let i = 0; i < unit.elements.length; i++) target.push(unit.elements[i]);
  }
  result.spheres.sort((a, b) => a - b);
  result.gaussians.sort((a, b) => a - b);
  return result;
}
function unitElements(structure: Structure) {
  return Object.fromEntries(structure.units.map((u) => [u.id, Array.from(u.elements)]));
}

const NSeq = 1337;
const Chains = ['A', 'B', 'C', 'D'];
const EntityOfChain: Record<string, string> = { A: '1', B: '1', C: '2', D: '2' };

type Category = {
  name: string;
  rows: MVSAnnotationRow[];
};

function makeCategories(rnd: () => number): Category[] {
  const chain = () => Chains[randInt(rnd, 4)];
  const seq = () => 1 + randInt(rnd, NSeq);
  const span = () => {
    const beg = seq();
    return { beg_label_seq_id: beg, end_label_seq_id: Math.min(NSeq, beg + randInt(rnd, 120)) };
  };
  const many = (n: number, make: () => MVSAnnotationRow) => Array.from({ length: n }, make);
  return [
    { name: 'chain only', rows: many(8, () => ({ label_asym_id: chain() })) },
    { name: 'entity only', rows: many(4, () => ({ label_entity_id: `${1 + randInt(rnd, 2)}` })) },
    { name: 'chain+seq', rows: many(40, () => ({ label_asym_id: chain(), label_seq_id: seq() })) },
    { name: 'chain+range', rows: many(40, () => ({ label_asym_id: chain(), ...span() })) },
    {
      name: 'chain+one-sided range',
      rows: many(20, () =>
        rnd() < 0.5
          ? { label_asym_id: chain(), beg_label_seq_id: seq() }
          : { label_asym_id: chain(), end_label_seq_id: seq() },
      ),
    },
    { name: 'seq only (all chains)', rows: many(15, () => ({ label_seq_id: seq() })) },
    { name: 'range only (all chains)', rows: many(15, () => span()) },
    { name: 'entity+range', rows: many(15, () => ({ label_entity_id: `${1 + randInt(rnd, 2)}`, ...span() })) },
    {
      name: 'entity+chain (contradicting about half of the time)',
      rows: many(20, () => {
        const c = chain();
        return { label_asym_id: c, label_entity_id: rnd() < 0.5 ? EntityOfChain[c] : `${1 + randInt(rnd, 2)}` };
      }),
    },
    {
      name: 'seq+range together (both have to hold)',
      rows: many(20, () => ({ label_asym_id: chain(), label_seq_id: seq(), ...span() })),
    },
    {
      name: 'inverted range (beg > end)',
      rows: [
        ...many(10, () => ({
          label_asym_id: chain(),
          beg_label_seq_id: 800 + randInt(rnd, 100),
          end_label_seq_id: 100 + randInt(rnd, 100),
        })),
        // (close bounds: an element which starts at 1 and ends at 20 contains both, but it is not a range)
        { label_asym_id: 'A', beg_label_seq_id: 20, end_label_seq_id: 1 },
        { label_asym_id: 'C', beg_label_seq_id: 2, end_label_seq_id: 1 },
        // (with a residue)
        { label_asym_id: 'A', beg_label_seq_id: 70, end_label_seq_id: 58, label_seq_id: 60 },
        { beg_label_seq_id: 30, end_label_seq_id: 10, label_seq_id: 20 },
        ...many(10, () => {
          const end = seq();
          return {
            label_asym_id: chain(),
            beg_label_seq_id: end + 1 + randInt(rnd, 30),
            end_label_seq_id: end,
            label_seq_id: seq(),
          };
        }),
      ],
    },
    {
      name: 'outside of the sequence',
      rows: [
        { label_asym_id: 'A', label_seq_id: NSeq + 5 },
        { label_asym_id: 'A', beg_label_seq_id: NSeq + 1 },
        { label_seq_id: 0 },
        { label_asym_id: 'Z' },
      ],
    },
    {
      name: 'fields which do not exist for coarse elements (select nothing)',
      rows: [
        { auth_asym_id: 'A' },
        { label_asym_id: 'A', auth_seq_id: 5 },
        { label_asym_id: 'B', label_atom_id: 'CA' },
        { label_comp_id: 'ALA' },
        { residue_index: 3 },
        { atom_id: 1 },
        { atom_index: 1 },
        { type_symbol: 'C' },
        { pdbx_PDB_ins_code: 'A' },
        { label_asym_id: 'C', beg_auth_seq_id: 1 },
        { auth_comp_id: 'GLY', label_asym_id: 'D' },
      ],
    },
    {
      name: 'instance_id',
      rows: [
        { instance_id: '1_555' },
        { instance_id: '1_555', label_asym_id: 'B' },
        { instance_id: 'other' },
        { instance_id: 'other', label_asym_id: 'A' },
        { operator_name: '1_555' },
        { operator_name: '1_555', label_asym_id: 'C' },
        { operator_name: 'other' },
        { operator_name: '1_555', instance_id: 'other' },
      ],
    },
    { name: 'no criteria', rows: [{}] },
  ];
}

describe('coarse models (spheres and gaussians)', async () => {
  const model = await loadCoarseModel();
  const structure = Structure.ofModel(model);
  const categories = makeCategories(mulberry32(2024));
  const oracle = createCoarseOracle(model);
  const allRows = categories.flatMap((c) => c.rows);

  it('the test structure is coarse only, with spheres and gaussians in several chains and entities', () => {
    expect(model.atomicHierarchy.atoms._rowCount).toBe(0);
    expect(model.coarseHierarchy.spheres.count).toBe(1764);
    expect(model.coarseHierarchy.gaussians.count).toBe(1764);
    expect(structure.units.map((u) => u.kind).sort()).toEqual(
      [Unit.Kind.Spheres, Unit.Kind.Spheres, Unit.Kind.Gaussians, Unit.Kind.Gaussians].sort(),
    );
  });

  const quiet = <T>(f: () => T): T => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    try {
      return f();
    } finally {
      warn.mockRestore();
    }
  };

  for (const c of categories) {
    it(`substructureFromRows selects what the definition of the row says: ${c.name}`, () => {
      quiet(() => expect(selectedElements(substructureFromRows(structure, c.rows))).toEqual(oracle.elementsOf(c.rows)));
    });
  }

  it('substructureFromRows for all categories together', () => {
    quiet(() => {
      const result = substructureFromRows(structure, allRows);
      expect(selectedElements(result)).toEqual(oracle.elementsOf(allRows));
      expect(result.elementCount).toBeGreaterThan(100); // (not vacuous)
    });
  });

  describe('getRowForEachAtom: the last row which selects an element wins', () => {
    const mapping = (rows: readonly MVSAnnotationRow[]) =>
      (new MVSAnnotation({ format: 'json', data: rows as any }, 'all_atomic', {}) as any).getRowForEachAtom(
        model,
        IdentityOperator,
      ) as {
        atoms: Int32Array | null;
        spheres: Int32Array | null;
        gaussians: Int32Array | null;
      };

    for (const c of categories) {
      it(c.name, () => {
        quiet(() => {
          const got = mapping(c.rows);
          const exp = oracle.rowOfEach(c.rows);
          expect(got.atoms).toBeNull();
          expect(Array.from(got.spheres ?? new Array(1764).fill(-1))).toEqual(exp.spheres);
          expect(Array.from(got.gaussians ?? new Array(1764).fill(-1))).toEqual(exp.gaussians);
        });
      });
    }
    it('all categories together', () => {
      quiet(() => {
        const got = mapping(allRows);
        const exp = oracle.rowOfEach(allRows);
        expect(Array.from(got.spheres!)).toEqual(exp.spheres);
        expect(Array.from(got.gaussians!)).toEqual(exp.gaussians);
      });
    });
  });

  describe('rows which cannot select a coarse element', () => {
    const inverted = categories.find((c) => c.name.startsWith('inverted'))!.rows;
    const atomicOnly = categories.find((c) => c.name.startsWith('fields which do not exist'))!.rows;
    const mapping = (rows: readonly MVSAnnotationRow[]) =>
      (new MVSAnnotation({ format: 'json', data: rows as any }, 'all_atomic', {}) as any).getRowForEachAtom(
        model,
        IdentityOperator,
      );

    it('a range with beg > end selects nothing, also when the bounds are close or a residue is given', () => {
      expect(substructureFromRows(structure, inverted).units.length).toBe(0);
      const got = mapping(inverted);
      expect(got.spheres).toBeNull();
      expect(got.gaussians).toBeNull();
      // (the sphere 1-20 contains the bounds 20 and 1, but is not in a range from 20 to 1)
      expect(
        substructureFromRows(structure, [{ label_asym_id: 'A', beg_label_seq_id: 20, end_label_seq_id: 1 }])
          .elementCount,
      ).toBe(0);
      expect(
        substructureFromRows(structure, [{ label_asym_id: 'A', beg_label_seq_id: 1, end_label_seq_id: 20 }])
          .elementCount,
      ).toBeGreaterThan(0);
    });
    it('values of another type, and a residue which is not a whole number, select nothing (as for atoms)', () => {
      const count = (row: object) => substructureFromRows(structure, [row as MVSAnnotationRow]).elementCount;
      expect(count({ label_asym_id: 'A', label_seq_id: 20 })).toBeGreaterThan(0);
      for (const row of [
        { label_asym_id: 'A', label_seq_id: '20' },
        { label_asym_id: 'A', beg_label_seq_id: '1' },
        { label_asym_id: 'A', end_label_seq_id: '100' },
        { label_asym_id: 1 },
        { label_entity_id: 1 },
        { label_asym_id: 'A', label_seq_id: 20.5 },
      ]) {
        expect(count(row), JSON.stringify(row)).toBe(0);
      }
      // (bounds need not be whole numbers)
      expect(count({ label_asym_id: 'A', beg_label_seq_id: 20.5, end_label_seq_id: 20.6 })).toBe(
        count({ label_asym_id: 'A', label_seq_id: 21 }),
      );
      for (const missing of [undefined, null, Number.NaN]) {
        expect(count({ label_asym_id: 'A', label_seq_id: missing })).toBe(count({ label_asym_id: 'A' }));
      }
    });
    it('criteria which only atoms have (auth_*, comp_id, atom ids, ...) select nothing, and warn once', () => {
      const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
      try {
        expect(substructureFromRows(structure, atomicOnly).units.length).toBe(0);
        const got = mapping(atomicOnly);
        expect(got.spheres).toBeNull();
        expect(got.gaussians).toBeNull();
        expect(substructureFromRows(structure, [{ auth_asym_id: 'A' }]).units.length).toBe(0); // (does not throw)
        expect(warn.mock.calls.length).toBeLessThanOrEqual(1); // (once per session; it may have been printed by an earlier test)
      } finally {
        warn.mockRestore();
      }
    });
  });
});
