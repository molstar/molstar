/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { SortedArray } from '@molstar/core/data/int';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { MoleculeType } from '@molstar/model/model/structure/model/types';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { computeInteractions } from '../interactions.js';
import { BridgeProviderParams, ContactProviderParams, type InteractionsProps } from '../params.js';

// Pins the inter-unit contacts (pairs, types, flags and their order) and bridges of symmetry-expanded structures, as
// computed before the inter-unit loop fetched A's features from A's grid (06-model-props-script N1). With every
// provider on, the flags changed once since: the refiners filter inter-unit contacts as intra-unit ones (the inter-unit
// check never matched), and a salt bridge hides the H-bonds within it rather than itself (as NGL and PLIP). Against the
// first pins below: 82, 4 and 54 inter-unit H-bonds more flagged; no ionic contact flagged in either.

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');

async function loadCif(name: string): Promise<Model> {
  const parsed = await CIF.parse(fs.readFileSync(path.join(examples, name), 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
}

/** Splits every unit at residue boundaries into chunks of about `size` elements, so contacts cross units. */
function splitUnits(s: Structure, size: number): Structure {
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
      const e = elements[i];
      if (chunk.length >= size && (!residueIndex || residueIndex[e] !== residueIndex[elements[i - 1]])) flush();
      chunk.push(e);
    }
    flush();
  }
  return Structure.create(units);
}

function cells(s: Structure, max: Vec3) {
  return StructureSymmetry.buildSymmetryRange(s, Vec3.create(0, 0, 0), max).run();
}

function hash(xs: ArrayLike<number>) {
  return (hashFnv32a(xs) >>> 0).toString(16);
}

function compareRows(a: number[], b: number[]) {
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return a[i] - b[i];
  return 0;
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

async function signature(s: Structure, props: Partial<InteractionsProps>) {
  const ctx = { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() };
  const { contacts, bridges } = await computeInteractions(ctx, s, props);
  const rows = contacts.edges.map((e) => [e.unitA, e.indexA, e.unitB, e.indexB, e.props.type, e.props.flag]);
  const types = new Map<number, number>();
  for (const r of rows) types.set(r[4], (types.get(r[4]) ?? 0) + 1);
  return [
    `count=${rows.length}`,
    `types=${[...types].sort((a, b) => a[0] - b[0]).map(([t, c]) => `${t}:${c}`)}`,
    `ordered=${hash(rows.flat())}`,
    `sorted=${hash([...rows].sort(compareRows).flat())}`,
    `bridges=${bridges.length}`,
  ].join(' ');
}

describe('inter-unit contacts', async () => {
  const m1cbs = await loadCif('1cbs_updated.cif');
  const m1bna = await loadCif('1bna_confal_pyramids.cif');
  const m7qpd = await loadCif('7qpd.fw2.cif');

  const cases: [name: string, structure: () => Promise<Structure>, defaults: string, all: string][] = [
    // crystal contacts between symmetry copies sharing features
    [
      '1cbs 2x2x1 cells',
      () => cells(Structure.ofModel(m1cbs), Vec3.create(1, 1, 0)),
      'count=5802 types=4:5802 ordered=851c5001 sorted=1e4921ad bridges=0',
      'count=6628 types=1:62,4:5802,6:430,8:334 ordered=4243b3b1 sorted=f7df69d1 bridges=1327',
    ],
    // DNA copies
    [
      '1bna cell',
      () => cells(Structure.ofModel(m1bna), Vec3.create(1, 0, 0)),
      'count=1498 types=4:1498 ordered=3f8d3109 sorted=87123ec1 bridges=0',
      'count=1604 types=4:1498,8:106 ordered=6ac3c81 sorted=f4b97265 bridges=33',
    ],
    // identity operators (the NAG N features gone, 06/B11: the same contacts, feature indices shifted)
    [
      '7qpd',
      async () => Structure.ofModel(m7qpd),
      'count=94 types=2:6,4:88 ordered=d26d8de5 sorted=26280639 bridges=0',
      'count=200 types=1:10,2:6,4:88,6:70,8:26 ordered=9242e4b9 sorted=968a8ad5 bridges=0',
    ],
    // contacts within chains split into units, under symmetry operators
    [
      '1cbs split cell',
      () => cells(splitUnits(Structure.ofModel(m1cbs), 300), Vec3.create(1, 0, 0)),
      'count=3690 types=2:16,4:3674 ordered=26b72ab5 sorted=f2211271 bridges=0',
      'count=5254 types=1:44,2:16,4:3674,6:958,8:562 ordered=9e8d46e5 sorted=ad43f549 bridges=635',
    ],
  ];

  it.each(cases)('%s', async (_name, structure, defaults, all) => {
    expect(await signature(await structure(), {})).toBe(defaults);
    expect(await signature(await structure(), AllProviders)).toBe(all);
  });
});

describe('water bridges', async () => {
  // two models of one entry, e.g. an ensemble shown as all models
  const mA = await loadCif('1cbs_updated.cif');
  const mB = await loadCif('1cbs_updated.cif');
  const isWater = (u: Unit) =>
    Unit.isAtomic(u) &&
    u.model.atomicHierarchy.derived.residue.moleculeType[u.residueIndex[u.elements[0]]] === MoleculeType.Water;
  const WaterBridgesOnly: Partial<InteractionsProps> = {
    providers: Object.fromEntries(
      Object.keys(ContactProviderParams).map((k) => [k, { name: 'off', params: {} }]),
    ) as InteractionsProps['providers'],
    bridges: allOn(BridgeProviderParams) as InteractionsProps['bridges'],
  };

  /** The non-water units of `a` and the water units of `b`, renumbered */
  function mixed(a: Model, b: Model) {
    const units: Unit[] = [];
    const add = (u: Unit) =>
      units.push(
        Unit.create(
          units.length,
          units.length,
          units.length,
          u.traits,
          u.kind,
          u.model,
          u.conformation.operator,
          u.elements,
        ),
      );
    for (const u of Structure.ofModel(a).units) if (!isWater(u)) add(u);
    for (const u of Structure.ofModel(b).units) if (isWater(u)) add(u);
    return Structure.create(units);
  }

  it('do not join the atoms of one model through the water of another (perf-audit-2 06/B7)', async () => {
    const ctx = { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() };
    const same = await computeInteractions(ctx, mixed(mA, mA), WaterBridgesOnly);
    expect(same.bridges.length).toBeGreaterThan(0);
    const across = await computeInteractions(ctx, mixed(mA, mB), WaterBridgesOnly);
    expect(across.bridges.length).toBe(0);
  });
});
