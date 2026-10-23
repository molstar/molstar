/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { parsePDB } from '@molstar/io/reader/pdb/parser';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { trajectoryFromPDB } from '@molstar/model/formats/structure/pdb';
import { Structure } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import { InteractionFlag, InteractionType } from '../common.js';
import { computeInteractions } from '../interactions.js';
import { BridgeProviderParams, ContactProviderParams, type InteractionsProps } from '../params.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');
const ctx = { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() };

async function loadCif(name: string): Promise<Model> {
  const parsed = await CIF.parse(fs.readFileSync(path.join(examples, name), 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
}

async function loadPdbqt(name: string): Promise<Model> {
  const parsed = await parsePDB(fs.readFileSync(path.join(examples, name), 'utf8'), name, 'pdbqt').run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromPDB(parsed.result).run()).representative;
}

/** Every contact and bridge provider on, with its default params. */
function allOn(params: Record<string, PD.Mapped<any>>) {
  const ret: Record<string, { name: 'on'; params: any }> = {};
  for (const [k, p] of Object.entries(params)) {
    ret[k] = { name: 'on', params: PD.getDefaultValues((p.map('on') as PD.Group<any>).params) };
  }
  return ret;
}

const AllProviders: Partial<InteractionsProps> = {
  providers: allOn(ContactProviderParams) as InteractionsProps['providers'],
  bridges: allOn(BridgeProviderParams) as InteractionsProps['bridges'],
};

/** The refiners that filter a contact of `applicable` type when a contact of `filters` type joins the same atoms */
const OverlapRefiners: [applicable: InteractionType[], filters: InteractionType[]][] = [
  // salt bridges are kept, the H-bonds within them filtered (as NGL and PLIP)
  [[InteractionType.HydrogenBond, InteractionType.WeakHydrogenBond], [InteractionType.Ionic]],
  [[InteractionType.Hydrophobic, InteractionType.CationPi], [InteractionType.PiStacking]],
  [[InteractionType.Ionic], [InteractionType.MetalCoordination]],
];

describe('interactions refinement', async () => {
  const m1cbs = await loadCif('1cbs_updated.cif');
  const mAce2 = await loadPdbqt('ace2.pdbqt');

  it('filters inter- and intra-unit contacts that share atoms with a contact of the filtering type', async () => {
    const overlaps = { inter: 0, intra: 0, saltBridgeHBonds: 0 };
    for (const s of [
      await StructureSymmetry.buildSymmetryRange(
        Structure.ofModel(m1cbs),
        Vec3.create(0, 0, 0),
        Vec3.create(1, 1, 0),
      ).run(),
      Structure.ofModel(mAce2),
    ]) {
      const { contacts, unitsContacts, unitsFeatures } = await computeInteractions(ctx, s, AllProviders);
      const members = (unitId: number, feature: number) => {
        const { offsets, members } = unitsFeatures.get(unitId);
        return Array.from(members.subarray(offsets[feature], offsets[feature + 1]));
      };
      const shares = (xs: number[], ys: number[]) => xs.some((x) => ys.includes(x));
      const check = (
        kind: 'inter' | 'intra',
        type: InteractionType,
        flag: number,
        overlapsWith: (filters: InteractionType[]) => boolean,
      ) => {
        for (const [applicable, filters] of OverlapRefiners) {
          if (!applicable.includes(type) || !overlapsWith(filters)) continue;
          overlaps[kind]++;
          if (filters.includes(InteractionType.Ionic)) overlaps.saltBridgeHBonds++;
          expect(flag).toBe(InteractionFlag.Filtered);
        }
      };

      for (const e of contacts.edges) {
        check('inter', e.props.type, e.props.flag, (filters) =>
          contacts.edges.some(
            (f) =>
              f.unitA === e.unitA &&
              f.unitB === e.unitB &&
              filters.includes(f.props.type) &&
              shares(members(f.unitA, f.indexA), members(e.unitA, e.indexA)) &&
              shares(members(f.unitB, f.indexB), members(e.unitB, e.indexB)),
          ),
        );
      }

      for (const unit of s.units) {
        const c = unitsContacts.get(unit.id);
        if (!c) continue;
        const { a, b, edgeProps } = c;
        for (let i = 0; i < c.edgeCount * 2; ++i) {
          check('intra', edgeProps.type[i], edgeProps.flag[i], (filters) => {
            for (let j = 0; j < c.edgeCount * 2; ++j) {
              if (!filters.includes(edgeProps.type[j])) continue;
              if (
                shares(members(unit.id, a[j]), members(unit.id, a[i])) &&
                shares(members(unit.id, b[j]), members(unit.id, b[i]))
              )
                return true;
            }
            return false;
          });
        }
      }
    }
    // the structures have such contacts, salt bridges with an H-bond among them, so the test checks something
    expect(overlaps.inter).toBeGreaterThan(0);
    expect(overlaps.intra).toBeGreaterThan(0);
    expect(overlaps.saltBridgeHBonds).toBeGreaterThan(0);
  });

  it('refines weak H-bonds without intra-unit contacts (`skipIntraContacts`, the interactions extension)', async () => {
    const s = await StructureSymmetry.buildSymmetryRange(
      Structure.ofModel(m1cbs),
      Vec3.create(0, 0, 0),
      Vec3.create(1, 1, 0),
    ).run();
    const { contacts, unitsContacts } = await computeInteractions(ctx, s, AllProviders, { skipIntraContacts: true });
    expect(unitsContacts.size).toBe(0);
    const types = new Set(contacts.edges.map((e) => e.props.type));
    expect(types.has(InteractionType.WeakHydrogenBond)).toBe(true);
    expect(types.has(InteractionType.HydrogenBond)).toBe(true);
  });
});
