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
import { parsePDB } from '@molstar/io/reader/pdb/parser';
import { trajectoryFromPDB } from '@molstar/model/formats/structure/pdb';
import { Structure, type StructureElement, Unit } from '@molstar/model/model/structure';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { eachIntraBondedAtom, typeSymbol } from '../../chemistry/util.js';
import { FeatureType, InteractionType } from '../common.js';
import { computeInteractions } from '../interactions.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');
const ctx = { runtime: SyncRuntimeContext, assetManager: new AssetManager() };

/** Donor–acceptor element pairs of the intra- and inter-unit H-bonds */
async function hydrogenBonds(structure: Structure) {
  const { contacts, unitsContacts, unitsFeatures } = await computeInteractions(ctx, structure, {});
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
  return ret;
}

describe('hydrogen bonds', async () => {
  // polar hydrogens (PDBQT)
  const parsed = await parsePDB(fs.readFileSync(path.join(examples, 'ace2.pdbqt'), 'utf8'), 'ace2', 'pdbqt').run();
  if (parsed.isError) throw new Error(parsed.message);
  const model = (await trajectoryFromPDB(parsed.result).run()).representative;
  const structure = Structure.ofModel(model);

  it("are rejected when the donor's only hydrogen is not in its unit (no fallback to the heavy atom)", async () => {
    const unit = structure.units[0] as Unit.Atomic;
    const { unitsContacts, unitsFeatures } = await computeInteractions(ctx, structure, {});
    const contacts = unitsContacts.get(unit.id);
    const { members, offsets, types } = unitsFeatures.get(unit.id);

    // an H-bond whose donor has exactly one hydrogen
    let donor = -1;
    let acceptor = -1;
    let hydrogen = -1;
    for (let i = 0; i < contacts.edgeCount * 2 && donor < 0; ++i) {
      if (contacts.edgeProps.type[i] !== InteractionType.HydrogenBond) continue;
      if (types[contacts.a[i]] !== FeatureType.HydrogenDonor) continue;
      const d = members[offsets[contacts.a[i]]] as StructureElement.UnitIndex;
      const hs: number[] = [];
      eachIntraBondedAtom(unit, d, (_, x) => {
        if (typeSymbol(unit, x) === 'H') hs.push(x);
      });
      if (hs.length !== 1) continue;
      donor = unit.elements[d];
      acceptor = unit.elements[members[offsets[contacts.b[i]]]];
      hydrogen = unit.elements[hs[0]];
    }
    expect(donor).toBeGreaterThanOrEqual(0);
    expect(await hydrogenBonds(structure)).toContain(`${donor}-${acceptor}`);

    // that hydrogen in a unit of its own: bonded to the donor across units
    const rest = SortedArray.ofSortedArray<ElementIndex>(
      Int32Array.from(unit.elements.filter((e) => e !== hydrogen)) as any,
    );
    const units: Unit[] = [Unit.create(0, 0, 0, unit.traits, unit.kind, model, unit.conformation.operator, rest)];
    const single = SortedArray.ofSingleton(hydrogen) as SortedArray<ElementIndex>;
    units.push(Unit.create(1, 1, 1, unit.traits, unit.kind, model, unit.conformation.operator, single));
    const split = Structure.create(units);

    const bonds = await hydrogenBonds(split);
    expect(bonds.size).toBeGreaterThan(0);
    expect(bonds).not.toContain(`${donor}-${acceptor}`);
  });
});
