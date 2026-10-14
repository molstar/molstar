/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { SymmetryOperator } from '@molstar/core/math/geometry/symmetry-operator';
import { Structure } from '@molstar/model/model/structure';
import type { MVSAnnotationRow } from '../schemas.js';
import { substructureFromRows } from '@molstar/model/model/structure/structure/element/schema-selection/row-structure';
import { combinedCategories, createModel, makeRows, rowCategories } from './annotation-rows-fixture.js';
import { createRowOracle } from './row-oracle.js';

/** What the structure looks like for comparisons: unit id -> its elements. */
function unitElements(structure: Structure) {
  return Object.fromEntries(structure.units.map((u) => [u.id, Array.from(u.elements)]));
}

/** The reference implementation: the atoms which satisfy a row, by the per-atom definition (`row-oracle.ts`), in the units of the structure. */
function reference(
  oracle: ReturnType<typeof createRowOracle>,
  structure: Structure,
  rows: readonly MVSAnnotationRow[],
) {
  const byOperator = new Map<string, Set<number>>();
  const result: Record<number, number[]> = {};
  for (const unit of structure.units) {
    const { instanceId, name } = unit.conformation.operator;
    const key = `${instanceId}|${name}`;
    let atoms = byOperator.get(key);
    if (!atoms) byOperator.set(key, (atoms = new Set(oracle.atomsOf(rows, { instanceId, name }))));
    const elements = Array.from(unit.elements).filter((e) => atoms.has(e));
    if (elements.length) result[unit.id] = elements;
  }
  return result;
}

/** Two copies of every unit of `structure`, the second one with a different operator (`instanceId` `copy-2`). */
function withSecondInstance(structure: Structure) {
  const builder = Structure.Builder({});
  const operator = SymmetryOperator.create('2_555', Mat4.fromTranslation(Mat4(), Vec3.create(100, 0, 0)), {
    instanceId: 'copy-2',
  });
  for (const unit of structure.units) builder.copyUnit(unit);
  for (const unit of structure.units) builder.addWithOperator(unit, operator);
  return builder.getStructure();
}

describe.each([
  [
    'single comp_id per residue',
    { nChains: 10, nResidues: 30, insertionCodes: true, descendingSeqIds: true, differentAuthComp: true },
  ],
  [
    'mixed comp_ids within residues',
    {
      nChains: 10,
      nResidues: 30,
      insertionCodes: true,
      descendingSeqIds: true,
      differentAuthComp: true,
      mixedCompIds: true,
    },
  ],
  [
    'repeated atom ids',
    {
      nChains: 12,
      nResidues: 30,
      insertionCodes: true,
      descendingSeqIds: true,
      differentAuthComp: true,
      atomIdWrap: 150,
    },
  ],
] as const)('substructureFromRows vs per-atom brute force (%s)', async (_label, dims) => {
  const model = await createModel(dims);
  const oracle = createRowOracle(model);
  const nAtoms = model.atomicHierarchy.atoms._rowCount;
  const whole = Structure.ofModel(model);
  const twoInstances = withSecondInstance(whole);
  const everyOtherUnit = (() => {
    const b = whole.subsetBuilder(true);
    whole.units.forEach((u, i) => {
      if (i % 2 === 0) b.setUnit(u.id, u.elements as any);
    });
    return b.getStructure();
  })();
  const instanceIds = Array.from(new Set(twoInstances.units.map((u) => u.conformation.operator.instanceId)));
  const operatorNames = Array.from(new Set(twoInstances.units.map((u) => u.conformation.operator.name)));

  const groups = [
    ...makeRows(rowCategories(dims, nAtoms, 400), 31),
    ...makeRows(combinedCategories(dims, nAtoms, 25), 32),
    // rows which pick one of the instances by their real `instance_id`
    {
      name: 'instance_id (real ids)',
      rows: instanceIds.flatMap((id) => [
        { instance_id: id },
        { instance_id: id, label_asym_id: 'B' },
        { instance_id: id, label_asym_id: 'C', label_seq_id: 3 },
      ]),
    },
    // ... and by the real name of the operator, alone or together with an instance (also one which belongs to the other operator)
    {
      name: 'operator_name (real names)',
      rows: operatorNames.flatMap((name) => [
        { operator_name: name },
        { operator_name: name, label_asym_id: 'B' },
        ...instanceIds.map((id) => ({ operator_name: name, instance_id: id, label_asym_id: 'D' })),
      ]),
    },
  ];

  it('the two-instance structure really has two operators', () => {
    expect(instanceIds.length).toBe(2);
    expect(operatorNames.sort()).toEqual(['1_555', '2_555']);
    expect(twoInstances.units.length).toBe(2 * whole.units.length);
  });

  const structures = [
    ['whole', whole],
    ['two instances', twoInstances],
    ['every other unit', everyOtherUnit],
  ] as const;

  for (const [structureName, structure] of structures) {
    for (const g of groups) {
      it(`${structureName}: ${g.name}`, () => {
        const result = substructureFromRows(structure, g.rows);
        const expected = reference(oracle, structure, g.rows);
        expect(unitElements(result)).toEqual(expected);
        if (result !== structure && result.elementCount > 0) expect(result.parent).toBe(structure.parent ?? structure);
      });
    }
    it(`${structureName}: all categories together`, () => {
      const rows = groups.flatMap((g) => g.rows);
      expect(unitElements(substructureFromRows(structure, rows))).toEqual(reference(oracle, structure, rows));
    });
  }

  it('no rows -> empty structure', () => {
    expect(substructureFromRows(whole, []).units.length).toBe(0);
    expect(substructureFromRows(Structure.Empty, [{ label_asym_id: 'A' }]).units.length).toBe(0);
  });

  it('selecting everything returns the source structure itself', () => {
    expect(substructureFromRows(whole, [{}])).toBe(whole);
    expect(substructureFromRows(everyOtherUnit, [{}])).toBe(everyOtherUnit);
  });
});
