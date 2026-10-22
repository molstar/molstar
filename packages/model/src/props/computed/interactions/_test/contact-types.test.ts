/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { SortedArray } from '@molstar/core/data/int';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, Unit } from '@molstar/model/model/structure';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import { CationPiProvider, IonicProvider, PiStackingProvider } from '../charged.js';
import type { FeatureType } from '../common.js';
import { addStructureContacts, addUnitContacts, type ContactProvider, type ContactTester } from '../contacts.js';
import { InterContactsBuilder, IntraContactsBuilder } from '../contacts-builder.js';
import type { Features } from '../features.js';
import { HalogenBondsProvider } from '../halogen-bonds.js';
import { HydrogenBondsProvider, WeakHydrogenBondsProvider } from '../hydrogen-bonds.js';
import { HydrophobicProvider } from '../hydrophobic.js';
import { computeInteractions } from '../interactions.js';
import { MetalCoordinationProvider } from '../metal.js';
import { BridgeProviderParams, ContactProviderParams, ContactsParams, type InteractionsProps } from '../params.js';

// The contact loops skip feature type pairs a tester does not accept (`ContactTester.acceptsTypes`) and query each
// feature only as far as its type's partners need (perf-audit-2 06/P4, 06/P10). With testers that accept every pair,
// the loops are what they were before; the contacts, in their order, must be the same.

const Providers: ContactProvider<any>[] = [
  IonicProvider,
  PiStackingProvider,
  CationPiProvider,
  HalogenBondsProvider,
  HydrogenBondsProvider,
  WeakHydrogenBondsProvider,
  HydrophobicProvider,
  MetalCoordinationProvider,
];

const testers = () => Providers.map((p) => p.createTester(PD.getDefaultValues(p.params)));

/** the same testers, accepting every type pair: the loops without the type masks */
const everyPair = (ts: ContactTester[]) => ts.map((t) => ({ ...t, acceptsTypes: () => true }));

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');

async function load1cbs() {
  const parsed = await CIF.parse(fs.readFileSync(path.join(examples, '1cbs_updated.cif'), 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  return Structure.ofModel((await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative);
}

/** Every unit split at residue boundaries into chunks of about 300 elements, under a symmetry operator or not. */
function split(s: Structure, symmetry: boolean) {
  const units: Unit[] = [];
  let id = 0;
  for (const u of s.units) {
    const { elements } = u;
    const residueIndex = Unit.isAtomic(u) ? u.residueIndex : undefined;
    let chunk: number[] = [];
    const flush = () => {
      if (chunk.length === 0) return;
      const set = SortedArray.ofSortedArray<ElementIndex>(Int32Array.from(chunk) as any);
      units.push(Unit.create(id, id, id, u.traits, u.kind, u.model, u.conformation.operator, set));
      id++;
      chunk = [];
    };
    for (let i = 0; i < elements.length; i++) {
      if (chunk.length >= 300 && (!residueIndex || residueIndex[elements[i]] !== residueIndex[elements[i - 1]])) {
        flush();
      }
      chunk.push(elements[i]);
    }
    flush();
  }
  const ret = Structure.create(units);
  // two of the space group's operators, contacts between them
  return symmetry ? StructureSymmetry.buildSymmetryRange(ret, Vec3.create(0, 0, 0), Vec3.create(0, 0, 0)).run() : ret;
}

function allOn(params: Record<string, PD.Mapped<any>>) {
  const ret: Record<string, { name: 'on'; params: any }> = {};
  for (const [k, p] of Object.entries(params)) {
    ret[k] = { name: 'on', params: PD.getDefaultValues((p.map('on') as PD.Group<any>).params) };
  }
  return ret;
}

/** The contacts of `s` for `ts` over the features every provider gives: intra rows per unit, then inter rows. */
function contactRows(s: Structure, unitsFeatures: { get(id: number): Features }, ts: ContactTester[]) {
  const props = PD.getDefaultValues(ContactsParams);
  const rows: number[][] = [];
  for (const u of s.units) {
    if (!Unit.isAtomic(u)) continue;
    const features = unitsFeatures.get(u.id);
    const builder = IntraContactsBuilder.create(features, u.elements.length);
    addUnitContacts(s, u, features, builder, ts, props);
    const c = builder.getContacts();
    for (let i = 0; i < c.edgeCount; i++) rows.push([u.id, c.a[i], c.b[i], c.edgeProps.type[i]]);
  }
  const builder = InterContactsBuilder.create();
  Structure.eachUnitPair(
    s,
    (a, b) => {
      const fa = unitsFeatures.get(a.id);
      const fb = unitsFeatures.get(b.id);
      addStructureContacts(s, a as Unit.Atomic, fa, b as Unit.Atomic, fb, builder, ts, props);
    },
    {
      maxRadius: Math.max(...ts.map((t) => t.maxDistance)),
      validUnit: (u) => Unit.isAtomic(u),
      validUnitPair: (a, b) => Structure.validUnitPair(s, a, b),
    },
  );
  for (const e of builder.getContacts(unitsFeatures as any).edges) {
    rows.push([e.unitA, e.indexA, e.unitB, e.indexB, e.props.type]);
  }
  return rows;
}

describe('contact testers by feature type pair', () => {
  it('give no type outside the pairs they accept', () => {
    for (const tester of testers()) {
      expect(tester.acceptsTypes).toBeDefined();
      for (let a = 0; a < 16; a++) {
        for (let b = 0; b < 16; b++) {
          if (tester.acceptsTypes!(a as FeatureType, b as FeatureType)) continue;
          // getType tests the types first: nothing else of the infos is read
          const info = (type: number) => ({ types: [type], feature: 0 }) as unknown as Features.Info;
          expect(tester.getType(undefined as any, info(a), info(b), 1)).toBeUndefined();
        }
      }
    }
  });

  it.each([false, true])('find the contacts of every pair, in the same order (symmetry %s)', async (symmetry) => {
    const s = await split(await load1cbs(), symmetry);
    const ctx = { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() };
    const props: Partial<InteractionsProps> = {
      providers: allOn(ContactProviderParams) as InteractionsProps['providers'],
      bridges: allOn(BridgeProviderParams) as InteractionsProps['bridges'],
    };
    const { unitsFeatures } = await computeInteractions(ctx, s, props);
    const ts = testers();
    const rows = contactRows(s, unitsFeatures, ts);
    expect(rows.length).toBeGreaterThan(500);
    expect(rows).toEqual(contactRows(s, unitsFeatures, everyPair(ts)));
  });
});
