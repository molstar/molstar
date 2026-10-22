/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { OrderedSet, SortedArray } from '@molstar/core/data/int';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, Unit } from '@molstar/model/model/structure';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import { CationPiProvider, IonicProvider, PiStackingProvider } from '../charged.js';
import { type FeatureType, InteractionType } from '../common.js';
import {
  addStructureContacts,
  addUnitContacts,
  checkLineOfSight,
  type ContactProvider,
  type ContactTester,
  validPair,
} from '../contacts.js';
import { InterContactsBuilder, IntraContactsBuilder } from '../contacts-builder.js';
import { Features } from '../features.js';
import { HalogenBondsProvider } from '../halogen-bonds.js';
import { HydrogenBondsProvider, WeakHydrogenBondsProvider } from '../hydrogen-bonds.js';
import { HydrophobicProvider } from '../hydrophobic.js';
import { computeInteractions } from '../interactions.js';
import { MetalCoordinationProvider } from '../metal.js';
import { BridgeProviderParams, ContactProviderParams, ContactsParams, type InteractionsProps } from '../params.js';

// The contact loops skip feature type pairs a tester does not accept (`ContactTester.acceptsTypes`) and query each
// feature only as far as its type's partners need (perf-audit-2 06/P4, 06/P10). The contacts, in their order, must be
// those of the loops before that change (`Reference` below), and those of testers that accept every pair.

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

async function parse(data: string | Uint8Array) {
  const parsed = await CIF.parse(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  return Structure.ofModel((await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative);
}

const load1cbs = () => parse(fs.readFileSync(path.join(examples, '1cbs_updated.cif'), 'utf8'));
/** zinc sites (metal coordination between units) and stacked bases */
const load1tup = () => parse(new Uint8Array(fs.readFileSync(path.join(examples, '1tup.bcif'))));

/**
 * Two halogen bonds, C–Cl···O=C on the C–Cl axis with 120° at the oxygen, 3.1 Å: one within chain A (its residues 2
 * and 3), one between chains A and B. No fixture of the repository has a halogen bond.
 */
function halogenBonds() {
  const atoms: [
    asym: string,
    entity: number,
    seq: number,
    comp: string,
    atom: string,
    symbol: string,
    xyz: number[],
  ][] = [
    ['A', 1, 1, 'CLA', 'C1', 'C', [0, 0, 0]],
    ['A', 1, 1, 'CLA', 'CL1', 'CL', [1.75, 0, 0]],
    ['A', 1, 2, 'CLB', 'C1', 'C', [0, 10, 0]],
    ['A', 1, 2, 'CLB', 'CL1', 'CL', [1.75, 10, 0]],
    ['A', 1, 3, 'ACC', 'O1', 'O', [4.85, 10, 0]],
    ['A', 1, 3, 'ACC', 'C1', 'C', [5.465, 11.065, 0]],
    ['B', 2, 1, 'ACC', 'O1', 'O', [4.85, 0, 0]],
    ['B', 2, 1, 'ACC', 'C1', 'C', [5.465, 1.065, 0]],
  ];
  const rows = atoms.map(
    ([asym, entity, seq, comp, atom, symbol, [x, y, z]], i) =>
      `HETATM ${i + 1} ${symbol} ${atom} ${comp} ${asym} ${entity} . ${x} ${y} ${z} ${seq} ${comp} ${asym} ${atom} 1`,
  );
  return parse(`data_X
loop_
_entity.id
_entity.type
1 non-polymer
2 non-polymer
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
${rows.join('\n')}
`);
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

type ContactsProps = PD.Values<typeof ContactsParams>;

/**
 * The contact loops before perf-audit-2 06/P4: every tester over its `requiredFeatures` subset within a unit; between
 * units every feature of A (near B's sphere) queried at the longest tester distance, every hit through `validPair` and
 * the testers in order. They took A's features from its grid near B (`getCandidates`), in ascending order; scanning all
 * of them gives the same candidates in the same order.
 */
const Reference = {
  addUnitContacts(
    structure: Structure,
    unit: Unit.Atomic,
    features: Features,
    builder: IntraContactsBuilder,
    testers: ReadonlyArray<ContactTester>,
    props: ContactsProps,
  ) {
    for (const tester of testers) {
      const { x, y, z } = features;
      const { lookup3d, indices: subsetIndices } = features.subset(tester.requiredFeatures);
      const infoA = Features.Info(structure, unit, features);
      const infoB = { ...infoA };
      for (let t = 0, tl = OrderedSet.size(subsetIndices); t < tl; ++t) {
        const i = OrderedSet.getAt(subsetIndices, t);
        const { count, indices, squaredDistances } = lookup3d.find(x[i], y[i], z[i], tester.maxDistance);
        if (count === 0) continue;
        infoA.feature = i;
        for (let r = 0; r < count; ++r) {
          const j = OrderedSet.getAt(subsetIndices, indices[r]);
          if (j <= i) continue;
          infoB.feature = j;
          if (!validPair(structure, infoA, infoB)) continue;
          const type = tester.getType(structure, infoA, infoB, squaredDistances[r]);
          if (type && checkLineOfSight(structure, infoA, infoB, props.lineOfSightDistFactor)) builder.add(i, j, type);
        }
      }
    }
  },
  addStructureContacts(
    structure: Structure,
    unitA: Unit.Atomic,
    featuresA: Features,
    unitB: Unit.Atomic,
    featuresB: Features,
    builder: InterContactsBuilder,
    testers: ReadonlyArray<ContactTester>,
    props: ContactsProps,
  ) {
    const { x: xA, y: yA, z: zA } = featuresA;
    const { lookup3d } = featuresB;
    const imageTransform = Mat4.mul(Mat4(), unitB.conformation.operator.inverse, unitA.conformation.operator.matrix);
    const isNotIdentity = !Mat4.isIdentity(imageTransform);
    const imageA = Vec3();
    const maxDistance = Math.max(...testers.map((t) => t.maxDistance));
    const { center, radius } = lookup3d.boundary.sphere;
    const testDistanceSq = (radius + maxDistance) * (radius + maxDistance);
    const infoA = Features.Info(structure, unitA, featuresA);
    const infoB = Features.Info(structure, unitB, featuresB);
    builder.startUnitPair(unitA, unitB);
    for (let i = 0 as Features.FeatureIndex; i < featuresA.count; i++) {
      Vec3.set(imageA, xA[i], yA[i], zA[i]);
      if (isNotIdentity) Vec3.transformMat4(imageA, imageA, imageTransform);
      if (Vec3.squaredDistance(imageA, center) > testDistanceSq) continue;
      const { indices, count, squaredDistances } = lookup3d.find(imageA[0], imageA[1], imageA[2], maxDistance);
      if (count === 0) continue;
      infoA.feature = i;
      for (let r = 0; r < count; ++r) {
        const j = indices[r];
        infoB.feature = j;
        if (!validPair(structure, infoA, infoB)) continue;
        const distanceSq = squaredDistances[r];
        for (const tester of testers) {
          if (distanceSq < tester.maxDistance * tester.maxDistance) {
            const type = tester.getType(structure, infoA, infoB, distanceSq);
            if (type && checkLineOfSight(structure, infoA, infoB, props.lineOfSightDistFactor)) {
              builder.add(i, j, type);
              break;
            }
          }
        }
      }
    }
    builder.finishUnitPair();
  },
};

const Current = { addUnitContacts, addStructureContacts };

/** The contacts of `s` for `ts` over the features every provider gives: intra rows per unit, then inter rows. */
function contactRows(
  s: Structure,
  unitsFeatures: { get(id: number): Features },
  ts: ContactTester[],
  loops: typeof Current = Current,
) {
  const props = PD.getDefaultValues(ContactsParams);
  const rows: number[][] = [];
  for (const u of s.units) {
    if (!Unit.isAtomic(u)) continue;
    const features = unitsFeatures.get(u.id);
    const builder = IntraContactsBuilder.create(features, u.elements.length);
    loops.addUnitContacts(s, u, features, builder, ts, props);
    const c = builder.getContacts();
    for (let i = 0; i < c.edgeCount; i++) rows.push([u.id, c.a[i], c.b[i], c.edgeProps.type[i]]);
  }
  const builder = InterContactsBuilder.create();
  Structure.eachUnitPair(
    s,
    (a, b) => {
      const fa = unitsFeatures.get(a.id);
      const fb = unitsFeatures.get(b.id);
      loops.addStructureContacts(s, a as Unit.Atomic, fa, b as Unit.Atomic, fb, builder, ts, props);
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

  /**
   * Intra and inter type counts (`type:count`; inter contacts listed under both units), after checking the rows against
   * the old loops and testers that accept every pair
   */
  async function checkedTypes(s: Structure) {
    const ctx = { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() };
    const props: Partial<InteractionsProps> = {
      providers: allOn(ContactProviderParams) as InteractionsProps['providers'],
      bridges: allOn(BridgeProviderParams) as InteractionsProps['bridges'],
    };
    const { unitsFeatures } = await computeInteractions(ctx, s, props);
    const ts = testers();
    const rows = contactRows(s, unitsFeatures, ts);
    expect(rows).toEqual(contactRows(s, unitsFeatures, ts, Reference));
    expect(rows).toEqual(contactRows(s, unitsFeatures, everyPair(ts)));
    const counts = (rs: number[][]) => {
      const m = new Map<number, number>();
      for (const r of rs) m.set(r[r.length - 1], (m.get(r[r.length - 1]) ?? 0) + 1);
      return [...m].sort((x, y) => x[0] - y[0]).map(([t, c]) => `${t}:${c}`);
    };
    return { intra: counts(rows.filter((r) => r.length === 4)), inter: counts(rows.filter((r) => r.length === 5)) };
  }

  it.each([false, true])('find the contacts of the old loops, in the same order (1cbs, symmetry %s)', async (sym) => {
    const { intra, inter } = await checkedTypes(await split(await load1cbs(), sym));
    expect(intra.length + inter.length).toBeGreaterThan(3);
  });

  it('find the metal coordination and pi-stacking of the old loops (1tup)', async () => {
    // types: 1 ionic, 2 cation-pi, 3 pi-stacking, 4 H-bond, 6 hydrophobic, 7 metal coordination, 8 weak H-bond
    expect(InteractionType.MetalCoordination).toBe(7);
    expect(await checkedTypes(await load1tup())).toEqual({
      intra: ['3:46', '4:883', '8:320'],
      inter: ['1:10', '2:2', '4:1232', '6:8', '7:24', '8:158'],
    });
  });

  it('find the halogen bonds of the old loops, within and between units', async () => {
    expect(InteractionType.HalogenBond).toBe(5);
    // an inter-unit contact is listed under both of its units
    expect(await checkedTypes(await halogenBonds())).toEqual({ intra: ['5:1'], inter: ['5:2'] });
  });

  it('test every pair for a tester without acceptsTypes, and follow a reused testers array', async () => {
    const s = await split(await load1cbs(), true);
    const ctx = { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() };
    const { unitsFeatures } = await computeInteractions(ctx, s, {});
    const hbonds = HydrogenBondsProvider.createTester(PD.getDefaultValues(HydrogenBondsProvider.params));
    const inter = (rows: number[][]) => rows.filter((r) => r.length === 5);
    const expected = inter(contactRows(s, unitsFeatures, [hbonds]));
    expect(expected.length).toBeGreaterThan(100);
    // the inter-unit loop gives such a tester every pair, whatever its `requiredFeatures` say
    const { acceptsTypes: _, ...custom } = hbonds;
    const donorsOnly = { ...custom, requiredFeatures: new Set([hbonds.requiredFeatures.values().next().value!]) };
    expect(inter(contactRows(s, unitsFeatures, [donorsOnly]))).toEqual(expected);
    // the same array with another tester
    const ts = [HydrophobicProvider.createTester(PD.getDefaultValues(HydrophobicProvider.params))];
    contactRows(s, unitsFeatures, ts);
    ts[0] = hbonds;
    expect(inter(contactRows(s, unitsFeatures, ts))).toEqual(expected);
    ts.push(donorsOnly);
    expect(() => contactRows(s, unitsFeatures, ts)).not.toThrow();
  });
});
