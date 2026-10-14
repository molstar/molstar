/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
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
  type: { name: 'schema', params: items },
  nullIfEmpty: true,
  label: '',
  ...extra,
});
const dispose = (b: any, cache: unknown) =>
  (StructureComponent as any).definition.dispose({ b, params: undefined, cache });

describe('StructureComponent with a schema', async () => {
  const structure = await createStructure();

  it('is the standard component transformer (so everything that works on components works on it)', () => {
    expect(StructureComponent.id).toBe('ms-plugin.structure-component');
  });

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
    expect(obj.description).toBeDefined();
    expect(obj.label).toBe('Component');
    const named = createStructureComponent(
      structure,
      params([{ label_asym_id: 'B' }], { label: 'Chain B' }),
      {} as any,
    ) as any;
    expect(named.label).toBe('Chain B');
    expect(named.data.elementCount).toBe(6 * 4);
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

  describe('update', () => {
    const U = StateTransformer.UpdateResult;
    const items = [{ label_asym_id: 'B' }];
    const make = () => {
      const cache: any = {};
      const obj = createStructureComponent(structure, params(items), cache) as any;
      return { cache, obj };
    };

    it('same schema, even as another object: Unchanged', () => {
      const { cache, obj } = make();
      expect(updateStructureComponent(structure, obj, params(items), params([{ label_asym_id: 'B' }]), cache)).toBe(
        U.Unchanged,
      );
    });
    it('only the label: Updated', () => {
      const { cache, obj } = make();
      expect(updateStructureComponent(structure, obj, params(items), params(items, { label: 'B' }), cache)).toBe(
        U.Updated,
      );
      expect(obj.label).toBe('B');
    });
    it('another schema, or another type: Recreate', () => {
      const { cache, obj } = make();
      expect(updateStructureComponent(structure, obj, params(items), params([{ label_asym_id: 'C' }]), cache)).toBe(
        U.Recreate,
      );
      expect(
        updateStructureComponent(
          structure,
          obj,
          params(items),
          { ...params(items), type: { name: 'static', params: 'all' } },
          cache,
        ),
      ).toBe(U.Recreate);
    });
    it('another structure with the same schema: the component is made from it, Updated (or Null if it became empty)', async () => {
      const { cache, obj } = make();
      const other = await createStructure();
      expect(updateStructureComponent(other, obj, params(items), params(items), cache)).toBe(U.Updated);
      expect(cache.source).toBe(other);
      expect(obj.data.elementCount).toBe(6 * 4);
      expect(obj.data.units.every((u: any) => other.unitMap.get(u.id) !== undefined)).toBe(true);
    });
    it('another structure which leaves nothing for the schema: Null', async () => {
      const cache: any = {};
      const obj = createStructureComponent(structure, params([{ label_asym_id: 'B' }]), cache) as any;
      const structureOfChainA = structure.subsetBuilder(true);
      for (const unit of structure.units) {
        if (unit.model.atomicHierarchy.chains.label_asym_id.value(unit.chainIndex[unit.elements[0]]) === 'A') {
          structureOfChainA.setUnit(unit.id, unit.elements as any);
        }
      }
      const withoutB = structureOfChainA.getStructure();
      expect(withoutB.elementCount).toBe(6 * 4);
      expect(
        updateStructureComponent(
          withoutB,
          obj,
          params([{ label_asym_id: 'B' }]),
          params([{ label_asym_id: 'B' }]),
          cache,
        ),
      ).toBe(U.Null);
    });
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
