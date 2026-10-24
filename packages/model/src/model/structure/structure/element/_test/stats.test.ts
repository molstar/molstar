/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Column, Table } from '@molstar/core/data/db';
import { OrderedSet, SortedArray } from '@molstar/core/data/int';
import { RuntimeContext } from '@molstar/core/task';
import { createModels } from '@molstar/model/formats/structure/basic/parser';
import { BasicSchema, createBasic } from '@molstar/model/formats/structure/basic/schema';
import { ComponentBuilder } from '@molstar/model/formats/structure/common/component';
import { EntityBuilder } from '@molstar/model/formats/structure/common/entity';
import { MoleculeType } from '../../../model/types.js';
import { Structure } from '../../structure.js';
import { Unit } from '../../unit.js';
import type { ElementIndex } from '../../../model/indexing.js';
import type { UnitIndex } from '../element.js';
import { Loci } from '../loci.js';
import { Stats } from '../stats.js';

/** Chains A (atoms 0, 1) and B (atoms 2, 3), one CA per residue */
async function twoChainModel() {
  const count = 4;
  const asym = Column.ofStringArray(['A', 'A', 'B', 'B']);
  const atom = Column.ofConst('CA', count, Column.Schema.str);
  const comp = Column.ofConst('ALA', count, Column.Schema.str);
  const seq = Column.ofIntArray([1, 2, 1, 2]);
  const atom_site = Table.ofPartialColumns(
    BasicSchema.atom_site,
    {
      auth_asym_id: asym,
      auth_atom_id: atom,
      auth_comp_id: comp,
      auth_seq_id: seq,
      Cartn_x: Column.ofFloatArray(Float32Array.of(0, 4, 8, 12)),
      Cartn_y: Column.ofFloatArray(new Float32Array(count)),
      Cartn_z: Column.ofFloatArray(new Float32Array(count)),
      id: Column.ofIntArray([1, 2, 3, 4]),
      label_asym_id: asym,
      label_atom_id: atom,
      label_comp_id: comp,
      label_seq_id: seq,
      label_entity_id: Column.ofConst('1', count, Column.Schema.str),
      occupancy: Column.ofConst(1, count, Column.Schema.float),
      type_symbol: Column.ofConst('C', count, Column.Schema.str),
      pdbx_PDB_model_num: Column.ofConst(1, count, Column.Schema.int),
    },
    count,
  );
  const entities = new EntityBuilder();
  entities.getEntityId('ALA', MoleculeType.Protein, 'A');
  const components = new ComponentBuilder(seq, atom);
  components.add('ALA', 0);
  const basic = createBasic({ entity: entities.getEntityTable(), chem_comp: components.getChemCompTable(), atom_site });
  const trajectory = await createModels(
    basic,
    { kind: 'test', name: 'two-chains', data: undefined },
    RuntimeContext.Synchronous,
  );
  return trajectory.representative;
}

describe('StructureElement.Stats', () => {
  it('finds the first whole chain across partitioned units (perf-audit-2 05/B4)', async () => {
    const model = await twoChainModel();
    const root = Structure.ofModel(model);
    const op = root.units[0].conformation.operator;
    const unit = (id: number, group: number, traits: Unit.Traits, elements: number[]) =>
      Unit.create(
        id,
        id,
        group,
        traits,
        Unit.Kind.Atomic,
        model,
        op,
        SortedArray.ofSortedArray(elements as ElementIndex[]),
      );
    // the first unit starts in chain A, which is not whole; chain B spans both units; the third keeps the loci from
    // being the whole structure
    const s = Structure.create([
      unit(0, 0, Unit.Trait.Partitioned | Unit.Trait.MultiChain, [1, 2]),
      unit(1, 0, Unit.Trait.Partitioned, [3]),
      unit(2, 1, Unit.Trait.None, [0]),
    ]);
    const loci = Loci(
      s,
      s.units
        .slice(0, 2)
        .map((u) => ({ unit: u, indices: OrderedSet.ofBounds(0 as UnitIndex, u.elements.length as UnitIndex) })),
    );
    const stats = Stats.ofLoci(loci);
    expect(stats.chainCount).toBe(1);
    expect(stats.firstChainLoc.element).toBe(2);
  });
});
