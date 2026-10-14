/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getHeapStatistics } from 'node:v8';
import { describe, expect, it, vi } from 'vitest';
import { SymmetryOperator } from '@molstar/core/math/geometry/symmetry-operator';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import { MVSAnnotation } from '../../components/annotation-prop.js';
import { substructureFromRows } from '@molstar/model/model/structure/structure/element/schema-selection/row-structure';
import type { MVSAnnotationRow } from '../schemas.js';
import { createModel } from './annotation-rows-fixture.js';
import { createRowOracle } from './row-oracle.js';

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

  const oracle = createRowOracle(model);
  /** The elements of each unit which satisfy the rows by the per-atom definition (units of the operators `op-a` and `op-b`, which have the same instance id) */
  const expectedUnits = (rows: MVSAnnotationRow[]) =>
    structure.units.flatMap((u) => {
      const { instanceId, name } = u.conformation.operator;
      const atoms = new Set(oracle.atomsOf(rows, { instanceId, name }));
      const elements = Array.from(u.elements).filter((e) => atoms.has(e));
      return elements.length ? [{ id: u.id, elements }] : [];
    });
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
      expect(result.units.map((u) => ({ id: u.id, elements: Array.from(u.elements) }))).toEqual(expectedUnits(rows));
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

  describe('the mapping of an annotation is computed once per what the rows depend on', () => {
    const locationIn = (unitIndex: number) => {
      const unit = structure.units[unitIndex];
      return StructureElement.Location.create(structure, unit, unit.elements[0]);
    };
    const indexOfA = structure.units.findIndex((u) => u.conformation.operator.name === 'op-a');
    const indexOfB = structure.units.findIndex((u) => u.conformation.operator.name === 'op-b');
    const make = (rows: any[]) => {
      const annot = new MVSAnnotation({ format: 'json', data: rows }, 'all_atomic', {});
      return { annot, compute: vi.spyOn(annot as any, 'getRowForEachAtom') };
    };

    it('rows which do not refer to an operator: once for all units', () => {
      const { annot, compute } = make([{ label_asym_id: 'A', tag: 'x' }]);
      for (let k = 0; k < 3; k++) {
        annot.getValueForLocation(locationIn(indexOfA), 'tag');
        annot.getValueForLocation(locationIn(indexOfB), 'tag');
      }
      expect(compute).toHaveBeenCalledTimes(1);
    });
    it('rows which refer to operator_name: once per operator, however many lookups', () => {
      const { annot, compute } = make([{ operator_name: 'op-a', label_asym_id: 'A', tag: 'x' }]);
      for (let k = 0; k < 5; k++) {
        annot.getValueForLocation(locationIn(indexOfA), 'tag');
        annot.getValueForLocation(locationIn(indexOfB), 'tag');
      }
      expect(compute).toHaveBeenCalledTimes(2);
    });
    it('rows which refer to an instance only: once per instance (these operators share theirs)', () => {
      const { annot, compute } = make([{ instance_id: 'shared-instance', label_asym_id: 'A', tag: 'x' }]);
      for (let k = 0; k < 3; k++) {
        annot.getValueForLocation(locationIn(indexOfA), 'tag');
        annot.getValueForLocation(locationIn(indexOfB), 'tag');
      }
      expect(compute).toHaveBeenCalledTimes(1);
    });
    it('the key of an operator is built when the operator is first looked up, not for every lookup', () => {
      for (const rows of [
        [{ instance_id: 'shared-instance', label_asym_id: 'A', tag: 'x' }],
        [{ operator_name: 'op-b', label_asym_id: 'A', tag: 'x' }],
      ]) {
        const { annot } = make(rows);
        const key = vi.spyOn(annot as any, 'operatorKey');
        for (let k = 0; k < 1000; k++) {
          annot.getValueForLocation(locationIn(indexOfA), 'tag');
          annot.getValueForLocation(locationIn(indexOfB), 'tag');
        }
        expect(key.mock.calls.length).toBeLessThanOrEqual(2); // (one per operator of the units)
      }
      const { annot } = make([{ label_asym_id: 'A', tag: 'x' }]);
      const key = vi.spyOn(annot as any, 'operatorKey');
      for (let k = 0; k < 100; k++) annot.getValueForLocation(locationIn(indexOfB), 'tag');
      expect(key).not.toHaveBeenCalled(); // (no row refers to an operator: there is no key)
    });
    // (`total_allocated_bytes` is in the heap statistics of newer versions of Node only)
    it.skipIf(typeof getHeapStatistics().total_allocated_bytes !== 'number')(
      'a lookup does not allocate (themes and tooltips look up every element)',
      () => {
        for (const rows of [
          [{ label_asym_id: 'A', tag: 'x' }],
          [{ instance_id: 'shared-instance', label_asym_id: 'A', tag: 'x' }],
          [{ operator_name: 'op-b', label_asym_id: 'A', tag: 'x' }],
        ]) {
          const { annot } = make(rows);
          const location = locationIn(indexOfB);
          for (let k = 0; k < 100000; k++) annot.getValueForLocation(location, 'tag'); // (builds the mapping, and warms the code up)
          const before = getHeapStatistics().total_allocated_bytes;
          for (let k = 0; k < 200000; k++) annot.getValueForLocation(location, 'tag');
          const perLookup = (getHeapStatistics().total_allocated_bytes - before) / 200000;
          expect(perLookup).toBeLessThan(8); // (a string built per lookup would be dozens of bytes)
        }
      },
    );
  });
});
