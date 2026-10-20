/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { SortedArray } from '@molstar/core/data/int';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { AssetManager } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { parsePDB } from '@molstar/io/reader/pdb/parser';
import { trajectoryFromPDB } from '@molstar/model/formats/structure/pdb';
import { Structure, Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { ValenceModelProvider } from '../../valence-model.js';
import { InteractionType } from '../common.js';
import { computeInteractions } from '../interactions.js';
import { ContactProviderParams, type InteractionsProps } from '../params.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');
const ctx = { runtime: SyncRuntimeContext, assetManager: new AssetManager() };

async function loadPdbqt(name: string): Promise<Model> {
  const parsed = await parsePDB(fs.readFileSync(path.join(examples, name), 'utf8'), name, 'pdbqt').run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromPDB(parsed.result).run()).representative;
}

describe('hydrogen bonds', async () => {
  // polar hydrogens (PDBQT)
  const mAce2 = await loadPdbqt('ace2.pdbqt');

  it('H-bond geometry follows a donor hydrogen in another unit', async () => {
    // every hydrogen of each unit moved to a unit of its own: donor N–H and O–H bonds cross units
    const s = Structure.ofModel(mAce2);
    const units: Unit[] = [];
    let id = 0;
    for (const u of s.units) {
      const { type_symbol } = u.model.atomicHierarchy.atoms;
      const heavy: number[] = [];
      const hydrogens: number[] = [];
      for (const e of u.elements) (type_symbol.value(e) === 'H' ? hydrogens : heavy).push(e);
      for (const set of [heavy, hydrogens]) {
        if (!set.length) continue;
        const elements = SortedArray.ofSortedArray<ElementIndex>(Int32Array.from(set) as any);
        units.push(Unit.create(id, id, id, u.traits, u.kind, u.model, u.conformation.operator, elements));
        id++;
      }
    }
    const split = Structure.create(units);
    expect(split.units.length).toBeGreaterThan(s.units.length);

    /** H-bonds by the elements of their features' first members, either unit layout */
    const hydrogenBonds = async (structure: Structure) => {
      // the valence model decides implicit hydrogens per unit (none if the unit has explicit ones), so a heavy-atom
      // unit would get donors the original does not have: keep to the explicit hydrogens
      await ValenceModelProvider.attach(ctx, structure, { assignH: 'never' });
      const props = { providers: { ...PD.getDefaultValues(ContactProviderParams) } } as Partial<InteractionsProps>;
      const { contacts, unitsContacts, unitsFeatures } = await computeInteractions(ctx, structure, props);
      const ret = new Set<string>();
      const element = (unit: Unit, feature: number) => {
        const { offsets, members } = unitsFeatures.get(unit.id);
        return unit.elements[members[offsets[feature]]];
      };
      for (const e of contacts.edges) {
        if (e.props.type !== InteractionType.HydrogenBond) continue;
        ret.add(
          `${element(structure.unitMap.get(e.unitA), e.indexA)}-${element(structure.unitMap.get(e.unitB), e.indexB)}`,
        );
      }
      for (const unit of structure.units) {
        const c = unitsContacts.get(unit.id);
        for (let i = 0; i < c.edgeCount * 2; ++i) {
          if (c.edgeProps.type[i] !== InteractionType.HydrogenBond) continue;
          ret.add(`${element(unit, c.a[i])}-${element(unit, c.b[i])}`);
        }
      }
      return [...ret].sort();
    };

    const expected = await hydrogenBonds(s);
    expect(expected.length).toBeGreaterThan(0);
    expect(await hydrogenBonds(split)).toEqual(expected);
  });
});
