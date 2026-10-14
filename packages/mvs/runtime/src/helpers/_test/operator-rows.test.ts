/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { SymmetryOperator } from '@molstar/core/math/geometry/symmetry-operator';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Structure, StructureElement, StructureSelection } from '@molstar/model/model/structure';
import { StructureQueryHelper } from '@molstar/plugin/state/helpers/structure-query';
import { MVSAnnotation } from '../../components/annotation-prop.js';
import { substructureFromRows } from '../row-structure.js';
import type { MVSAnnotationRow } from '../schemas.js';
import { rowsToExpression } from '../selections.js';
import { createModel } from './annotation-rows-fixture.js';

/** Two copies of `structure` whose operators have the same `instanceId` but different names, which rows can tell apart only by `operator_name`. */
function sameInstanceIdDifferentNames(structure: Structure) {
  const builder = Structure.Builder({});
  const shift = (name: string, x: number) =>
    SymmetryOperator.create(name, Mat4.fromTranslation(Mat4(), Vec3.create(x, 0, 0)), {
      instanceId: 'shared-instance',
    });
  for (const unit of structure.units) builder.addWithOperator(unit, shift('op-a', 100));
  for (const unit of structure.units) builder.addWithOperator(unit, shift('op-b', 200));
  return builder.getStructure();
}

describe('operator_name when the instance id does not tell the operators apart', async () => {
  const model = await createModel({ nChains: 3, nResidues: 10 });
  const structure = sameInstanceIdDifferentNames(Structure.ofModel(model));
  const nUnits = structure.units.length / 2;
  const unitNames = structure.units.map((u) => u.conformation.operator.name);

  it('the structure has the units of both operators under one instance id', () => {
    expect(new Set(structure.units.map((u) => u.conformation.operator.instanceId))).toEqual(
      new Set(['shared-instance']),
    );
    expect(unitNames.filter((n) => n === 'op-a').length).toBe(nUnits);
    expect(unitNames.filter((n) => n === 'op-b').length).toBe(nUnits);
  });

  const viaExpression = (rows: MVSAnnotationRow[]) =>
    StructureSelection.unionStructure(StructureQueryHelper.createAndRun(structure, rowsToExpression(rows)).selection);
  const unitNamesOf = (s: Structure) => Array.from(new Set(s.units.map((u) => u.conformation.operator.name))).sort();

  it('substructureFromRows selects the units of the named operator only', () => {
    for (const rows of [
      [{ operator_name: 'op-a' }],
      [{ operator_name: 'op-b', label_asym_id: 'B' }],
      [{ operator_name: 'op-a' }, { operator_name: 'op-b', label_asym_id: 'C', label_seq_id: 2 }],
      [{ instance_id: 'shared-instance', operator_name: 'op-b' }],
      [{ instance_id: 'shared-instance' }],
      [{ operator_name: 'nope' }],
    ] as MVSAnnotationRow[][]) {
      const result = substructureFromRows(structure, rows);
      const reference = viaExpression(rows);
      expect(result.elementCount).toBe(reference.elementCount);
      expect(result.units.map((u) => u.id)).toEqual(reference.units.map((u) => u.id));
    }
    expect(unitNamesOf(substructureFromRows(structure, [{ operator_name: 'op-a' }]))).toEqual(['op-a']);
    expect(
      unitNamesOf(substructureFromRows(structure, [{ operator_name: 'op-a' }, { operator_name: 'op-b' }])),
    ).toEqual(['op-a', 'op-b']);
  });

  it('an annotation maps the elements of each operator separately', () => {
    const annot = new MVSAnnotation(
      { format: 'json', data: [{ operator_name: 'op-b', label_asym_id: 'A', tag: 'b' }] as any },
      'all_atomic',
      {},
    );
    const valueIn = (unitIndex: number) => {
      const unit = structure.units[unitIndex];
      return annot.getValueForLocation(StructureElement.Location.create(structure, unit, unit.elements[0]), 'tag');
    };
    const indexOfA = structure.units.findIndex((u) => u.conformation.operator.name === 'op-a');
    const indexOfB = structure.units.findIndex((u) => u.conformation.operator.name === 'op-b');
    // (chain A is the first unit of each copy). The first lookup must not decide for the other operator.
    expect(valueIn(indexOfA)).toBeUndefined();
    expect(valueIn(indexOfB)).toBe('b');
    expect(valueIn(indexOfA)).toBeUndefined();
  });
});
