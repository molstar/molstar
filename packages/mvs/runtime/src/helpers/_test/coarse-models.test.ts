/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import * as fs from 'node:fs';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, StructureSelection, Unit } from '@molstar/model/model/structure';
import { StructureQueryHelper } from '@molstar/plugin/state/helpers/structure-query';
import { MVSAnnotation } from '../../components/annotation-prop.js';
import { substructureFromRows } from '@molstar/model/model/structure/structure/element/schema-selection/row-structure';
import type { MVSAnnotationRow } from '../schemas.js';
import { rowsToExpression } from '../selections.js';
import { IdentityOperator, mulberry32, randInt } from './annotation-rows-fixture.js';

/** `data/examples/ihm-coarse-8zzh.cif`: 4 coarse chains of 882 elements (A, C: spheres; B, D: gaussians), 2 entities, no atoms. */
async function loadCoarseModel() {
  const parsed = await CIF.parseText(fs.readFileSync('data/examples/ihm-coarse-8zzh.cif', 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return trajectory.representative;
}

function viaExpression(structure: Structure, rows: readonly MVSAnnotationRow[]) {
  const { selection } = StructureQueryHelper.createAndRun(structure, rowsToExpression(rows));
  return StructureSelection.unionStructure(selection);
}
function unitElements(structure: Structure) {
  return Object.fromEntries(structure.units.map((u) => [u.id, Array.from(u.elements)]));
}

const NSeq = 1337;
const Chains = ['A', 'B', 'C', 'D'];
const EntityOfChain: Record<string, string> = { A: '1', B: '1', C: '2', D: '2' };

type Category = {
  name: string;
  rows: MVSAnnotationRow[] /** Not comparable with MolScript, see the tests of `differsFromMolScript` */;
  differsFromMolScript?: boolean;
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
      name: 'seq+range together (as MolScript does)',
      rows: many(20, () => ({ label_asym_id: chain(), label_seq_id: seq(), ...span() })),
    },
    {
      name: 'inverted range (beg > end)',
      differsFromMolScript: true,
      rows: [
        ...many(10, () => ({
          label_asym_id: chain(),
          beg_label_seq_id: 800 + randInt(rnd, 100),
          end_label_seq_id: 100 + randInt(rnd, 100),
        })),
        // (close bounds: an element which starts at 1 and ends at 20 contains both, MolScript selects it)
        { label_asym_id: 'A', beg_label_seq_id: 20, end_label_seq_id: 1 },
        { label_asym_id: 'C', beg_label_seq_id: 2, end_label_seq_id: 1 },
        // (with a residue, which MolScript's overlap test then takes as a part of the elements which it selects)
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
      differsFromMolScript: true,
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
  const comparable = categories.filter((c) => !c.differsFromMolScript);
  const allRows = categories.flatMap((c) => c.rows);
  const comparableRows = comparable.flatMap((c) => c.rows);

  it('the test structure is coarse only, with spheres and gaussians in several chains and entities', () => {
    expect(model.atomicHierarchy.atoms._rowCount).toBe(0);
    expect(model.coarseHierarchy.spheres.count).toBe(1764);
    expect(model.coarseHierarchy.gaussians.count).toBe(1764);
    expect(structure.units.map((u) => u.kind).sort()).toEqual(
      [Unit.Kind.Spheres, Unit.Kind.Spheres, Unit.Kind.Gaussians, Unit.Kind.Gaussians].sort(),
    );
  });

  for (const c of comparable) {
    it(`substructureFromRows equals the MolScript expression: ${c.name}`, () => {
      const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
      try {
        const result = substructureFromRows(structure, c.rows);
        expect(unitElements(result)).toEqual(unitElements(viaExpression(structure, c.rows)));
      } finally {
        warn.mockRestore();
      }
    });
  }

  it('substructureFromRows equals the MolScript expression for all categories together', () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    try {
      const result = substructureFromRows(structure, comparableRows);
      expect(unitElements(result)).toEqual(unitElements(viaExpression(structure, comparableRows)));
      expect(result.elementCount).toBeGreaterThan(100); // (not vacuous)
      // the other rows select nothing, so adding them changes nothing
      expect(unitElements(substructureFromRows(structure, allRows))).toEqual(unitElements(result));
    } finally {
      warn.mockRestore();
    }
  });

  describe('getRowForEachAtom: the last row which selects an element wins', () => {
    /** Brute force: run the MolScript expression of each row separately and let later rows overwrite earlier ones. */
    function expected(rows: readonly MVSAnnotationRow[]) {
      const spheres = new Array<number>(model.coarseHierarchy.spheres.count).fill(-1);
      const gaussians = new Array<number>(model.coarseHierarchy.gaussians.count).fill(-1);
      rows.forEach((row, i) => {
        for (const unit of viaExpression(structure, [row]).units) {
          const target = unit.kind === Unit.Kind.Spheres ? spheres : gaussians;
          for (let k = 0; k < unit.elements.length; k++) target[unit.elements[k]] = i;
        }
      });
      return { spheres, gaussians };
    }
    const mapping = (rows: readonly MVSAnnotationRow[]) =>
      (new MVSAnnotation({ format: 'json', data: rows as any }, 'all_atomic', {}) as any).getRowForEachAtom(
        model,
        IdentityOperator,
      ) as {
        atoms: Int32Array | null;
        spheres: Int32Array | null;
        gaussians: Int32Array | null;
      };

    for (const c of comparable) {
      it(c.name, () => {
        const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
        try {
          const got = mapping(c.rows);
          const exp = expected(c.rows);
          expect(got.atoms).toBeNull();
          expect(Array.from(got.spheres ?? new Array(1764).fill(-1))).toEqual(exp.spheres);
          expect(Array.from(got.gaussians ?? new Array(1764).fill(-1))).toEqual(exp.gaussians);
        } finally {
          warn.mockRestore();
        }
      });
    }
    it('all categories together', () => {
      const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
      try {
        const got = mapping(comparableRows);
        const exp = expected(comparableRows);
        expect(Array.from(got.spheres!)).toEqual(exp.spheres);
        expect(Array.from(got.gaussians!)).toEqual(exp.gaussians);
      } finally {
        warn.mockRestore();
      }
    });
  });

  describe('where the selection differs from the MolScript expression on purpose', () => {
    const inverted = categories.find((c) => c.name.startsWith('inverted'))!.rows;
    const atomicOnly = categories.find((c) => c.name.startsWith('fields which do not exist'))!.rows;
    const mapping = (rows: readonly MVSAnnotationRow[]) =>
      (new MVSAnnotation({ format: 'json', data: rows as any }, 'all_atomic', {}) as any).getRowForEachAtom(
        model,
        '1_555',
      );

    it('a range with beg > end selects nothing, also when the bounds are close or a residue is given (MolScript selects elements which span them)', () => {
      expect(substructureFromRows(structure, inverted).units.length).toBe(0);
      const got = mapping(inverted);
      expect(got.spheres).toBeNull();
      expect(got.gaussians).toBeNull();
      // (and this is a difference: MolScript selects the sphere 1-20 for the bounds 20 and 1)
      expect(
        viaExpression(structure, [{ label_asym_id: 'A', beg_label_seq_id: 20, end_label_seq_id: 1 }]).elementCount,
      ).toBeGreaterThan(0);
    });
    it('criteria which only exist for atoms (auth_*, comp_id, atom ids, ...) select nothing, MolScript throws', () => {
      const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
      try {
        expect(substructureFromRows(structure, atomicOnly).units.length).toBe(0);
        const got = mapping(atomicOnly);
        expect(got.spheres).toBeNull();
        expect(got.gaussians).toBeNull();
      } finally {
        warn.mockRestore();
      }
      expect(() => viaExpression(structure, [{ auth_asym_id: 'A' }])).toThrow(
        'Property only available for atomic models',
      );
    });
  });
});
