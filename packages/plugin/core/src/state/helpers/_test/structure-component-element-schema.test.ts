/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import '@molstar/model/script/runtime/query/table'; // (registers the symbols of the queries)
import { CIF } from '@molstar/io/reader/cif';
import { StateObject, StateTransformer } from '@molstar/core/state';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import { StructureComponent } from '../../transforms/structure/selection.js';
import {
  type StructureComponentParams,
  createStructureComponent,
  updateStructureComponent,
} from '../structure-component.js';

/** Chains A, B, C of 6 residues with 4 atoms */
function mmcif() {
  const lines = [
    'data_T',
    'loop_',
    ...'group_PDB id type_symbol label_atom_id label_alt_id label_comp_id label_asym_id label_entity_id label_seq_id pdbx_PDB_ins_code Cartn_x Cartn_y Cartn_z occupancy B_iso_or_equiv auth_seq_id auth_comp_id auth_asym_id auth_atom_id pdbx_PDB_model_num'
      .split(' ')
      .map((c) => `_atom_site.${c}`),
  ];
  let id = 1;
  for (const asym of ['A', 'B', 'C']) {
    for (let r = 1; r <= 6; r++) {
      for (const name of ['N', 'CA', 'C', 'O']) {
        lines.push(`ATOM ${id++} ${name[0]} ${name} . ALA ${asym} 1 ${r} ? 0 0 0 1 10 ${r} ALA ${asym} ${name} 1`);
      }
    }
  }
  return `${lines.join('\n')}\n`;
}

async function createStructure() {
  const parsed = await CIF.parseText(mmcif()).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

const params = (
  items: StructureElement.SchemaItem[],
  extra: Partial<StructureComponentParams> = {},
): StructureComponentParams => ({
  type: { name: 'expression', params: StructureElement.Schema.toExpression({ items }) },
  nullIfEmpty: true,
  label: '',
  ...extra,
});
const dispose = (b: any, cache: unknown) =>
  (StructureComponent as any).definition.dispose({ b, params: undefined, cache });

/** The standard `StructureComponent` with the expression of a schema (the generator `element-schema`) */
describe('StructureComponent with the expression of a schema', async () => {
  const structure = await createStructure();

  it('contains the elements which satisfy the items', () => {
    const cache: any = {};
    const obj = createStructureComponent(
      structure,
      params([
        { label_asym_id: 'A', beg_label_seq_id: 2, end_label_seq_id: 3 },
        { label_asym_id: 'C', label_seq_id: 6 },
      ]),
      cache,
    ) as any;
    expect(obj.data.elementCount).toBe(2 * 4 + 4);
    expect(cache.source).toBe(structure);
    const named = createStructureComponent(
      structure,
      params([{ label_asym_id: 'B' }], { label: 'Chain B' }),
      {} as any,
    ) as any;
    expect(named.label).toBe('Chain B');
    expect(named.data.elementCount).toBe(6 * 4);
  });
  it('many rows (given by columns) are the same component as the elements of the rows', () => {
    // 12 rows: residues 1-3 of every chain and a residue of chain C, for every row a few atoms
    const rows = [
      ...['A', 'B', 'C'].flatMap((asym) => [1, 2, 3].map((label_seq_id) => ({ label_asym_id: asym, label_seq_id }))),
      { label_asym_id: 'C', label_seq_id: 6, label_atom_id: 'CA' },
      { label_asym_id: 'C', label_seq_id: 5, label_atom_id: 'N' },
      { label_asym_id: 'Z', label_seq_id: 1 },
    ];
    const p = params(rows);
    const expression = p.type.params as any;
    expect(expression.args[0].args.label_seq_id).toEqual([1, 2, 3, 1, 2, 3, 1, 2, 3, 6, 5, 1]); // (columns)
    const obj = createStructureComponent(structure, p, {} as any) as any;
    expect(obj.data.elementCount).toBe(9 * 4 + 1 + 1);
  });
  it('is null if empty (and nullIfEmpty), otherwise an empty structure', () => {
    expect(createStructureComponent(structure, params([{ label_asym_id: 'Z' }]), {} as any)).toBe(StateObject.Null);
    expect(createStructureComponent(structure, params([]), {} as any)).toBe(StateObject.Null);
    const empty = createStructureComponent(
      structure,
      params([{ label_asym_id: 'Z' }], { nullIfEmpty: false }),
      {} as any,
    ) as any;
    expect(empty.data.elementCount).toBe(0);
  });

  describe('dispose', () => {
    it('disposes the properties of a component which is a part of the structure', () => {
      const cache: any = {};
      const obj = createStructureComponent(structure, params([{ label_asym_id: 'B' }]), cache) as any;
      expect(obj.data.parent).toBe(structure.root);
      const spy = vi.spyOn(obj.data.customPropertyDescriptors, 'dispose');
      dispose(obj, cache);
      expect(spy).toHaveBeenCalledTimes(1);
    });
    it('disposes the properties of a component which is not the source (a copy of it), and not those of the source', () => {
      const cache: any = {};
      const all: StructureComponentParams = { type: { name: 'static', params: 'all' }, nullIfEmpty: true, label: '' };
      const obj = createStructureComponent(structure, all, cache) as any;
      expect(obj.data).toBe(structure);
      expect(cache.isSource).toBe(true);
      // (what the component is after its model was replaced: a copy of the source, which is its own)
      obj.data = obj.data.remapModel(structure.model);
      cache.isSource = false;
      const own = vi.spyOn(obj.data.customPropertyDescriptors, 'dispose');
      const theirs = vi.spyOn(structure.customPropertyDescriptors, 'dispose');
      dispose(obj, cache);
      expect(own).toHaveBeenCalledTimes(1);
      expect(theirs).not.toHaveBeenCalled();
    });
    it('does not dispose the shared empty structure', () => {
      const cache: any = {};
      const obj = createStructureComponent(
        structure,
        params([{ label_asym_id: 'Z' }], { nullIfEmpty: false }),
        cache,
      ) as any;
      expect(obj.data).toBe(Structure.Empty);
      const spy = vi.spyOn(Structure.Empty.customPropertyDescriptors, 'dispose');
      dispose(obj, cache);
      expect(spy).not.toHaveBeenCalled();
    });
    it('follows the source when an expression component becomes the source, or stops being it', async () => {
      const cache: any = {};
      const everything = params([{}]);
      const obj = createStructureComponent(structure, everything, cache) as any;
      expect(cache.isSource).toBe(true);
      const half = structure.subsetBuilder(true);
      for (const unit of structure.units) half.setUnit(unit.id, unit.elements.slice(0, 12) as any);
      const part = half.getStructure();
      expect(updateStructureComponent(part, obj, everything, everything, cache)).toBe(
        StateTransformer.UpdateResult.Updated,
      );
      expect(obj.data).toBe(part); // (everything of the new source: still the source itself)
      expect(cache.isSource).toBe(true);
    });
    it('does not dispose those of the structure when every element satisfies the schema (the component is the structure)', () => {
      const cache: any = {};
      const obj = createStructureComponent(structure, params([{}]), cache) as any;
      expect(obj.data).toBe(structure);
      const spy = vi.spyOn(structure.customPropertyDescriptors, 'dispose');
      dispose(obj, cache);
      expect(spy).not.toHaveBeenCalled();
    });
    it('does not dispose those of the source when the source is a part of another structure and every element satisfies the schema', () => {
      const half = structure.subsetBuilder(true);
      for (const unit of structure.units) half.setUnit(unit.id, unit.elements.slice(0, 12) as any);
      const source = half.getStructure(); // a part of `structure`: it has a parent
      expect(source.parent).toBe(structure);
      const cache: any = {};
      const obj = createStructureComponent(source, params([{}]), cache) as any;
      expect(obj.data).toBe(source);
      const spy = vi.spyOn(source.customPropertyDescriptors, 'dispose');
      dispose(obj, cache);
      expect(spy).not.toHaveBeenCalled();
    });
  });
});
