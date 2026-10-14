/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { SymmetryOperator } from '@molstar/core/math/geometry/symmetry-operator';
import { Structure, StructureSelection } from '@molstar/model/model/structure';
import { StructureQueryHelper } from '@molstar/plugin/state/helpers/structure-query';
import type { MVSAnnotationRow } from '../schemas.js';
import { substructureFromRows } from '@molstar/model/model/structure/structure/element/schema-selection/row-structure';
import { rowsToExpression } from '../selections.js';
import { combinedCategories, createModel, makeRows, rowCategories } from './annotation-rows-fixture.js';

/** What the structure looks like for comparisons: unit id -> its elements. */
function unitElements(structure: Structure) {
  return Object.fromEntries(structure.units.map((u) => [u.id, Array.from(u.elements)]));
}

/** The reference implementation: run the MolScript expression for the rows on every atom. */
function viaExpression(structure: Structure, rows: readonly MVSAnnotationRow[]) {
  const { selection } = StructureQueryHelper.createAndRun(structure, rowsToExpression(rows));
  return StructureSelection.unionStructure(selection);
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
] as const)('substructureFromRows vs MolScript expression (%s)', async (_label, dims) => {
  const model = await createModel(dims);
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
        const reference = viaExpression(structure, g.rows);
        expect(unitElements(result)).toEqual(unitElements(reference));
        if (result === structure) {
          // everything is selected (the expression path then returns an equal child structure, or the source itself for a row without criteria)
          expect(reference.elementCount).toBe(structure.elementCount);
        } else {
          expect(result.parent === undefined).toBe(reference.parent === undefined);
          if (result.parent) expect(result.parent).toBe(reference.parent);
        }
      });
    }
    it(`${structureName}: all categories together`, () => {
      const rows = groups.flatMap((g) => g.rows);
      expect(unitElements(substructureFromRows(structure, rows))).toEqual(unitElements(viaExpression(structure, rows)));
    });
  }

  it('no rows -> empty structure', () => {
    expect(substructureFromRows(whole, []).units.length).toBe(0);
    expect(substructureFromRows(Structure.Empty, [{ label_asym_id: 'A' }]).units.length).toBe(0);
  });

  it('selecting everything returns the source structure itself, like the expression path', () => {
    const rows = [{}]; // everything
    const result = substructureFromRows(whole, rows);
    const reference = viaExpression(whole, rows);
    expect(result.elementCount).toBe(reference.elementCount);
    expect(result).toBe(reference);
    expect(result).toBe(whole);
  });
});
