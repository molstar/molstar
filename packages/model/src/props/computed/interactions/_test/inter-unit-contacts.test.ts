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
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { SymmetryOperator } from '@molstar/core/math/geometry/symmetry-operator';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { AssetManager } from '@molstar/core/util/assets';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { computeInteractions, Interactions, type InterUnitContactsStats } from '../interactions.js';
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

/**
 * C_n of `structure` about the z axis through its centre shifted by `offset` along x (a closed point group: every pair
 * of neighbours is related by the same rotation, computed by a different chain of rounding); `extra` adds units
 */
function cyclic(
  structure: Structure,
  n: number,
  offset: number,
  extra?: (builder: Structure.StructureBuilder) => void,
) {
  const axis = Vec3.add(Vec3(), structure.boundary.sphere.center, Vec3.create(offset, 0, 0));
  const builder = Structure.Builder({ label: `C${n}` });
  for (let k = 0; k < n; k++) {
    const m = Mat4.fromTranslation(Mat4(), axis);
    Mat4.mul(m, m, Mat4.fromRotation(Mat4(), (2 * Math.PI * k) / n, Vec3.unitZ));
    Mat4.mul(m, m, Mat4.fromTranslation(Mat4(), Vec3.negate(Vec3(), axis)));
    const op = SymmetryOperator.create(`C${n}-${k}`, m);
    for (const unit of structure.units) builder.addWithOperator(unit, op);
  }
  extra?.(builder);
  return builder.getStructure();
}

function isWater(u: Unit) {
  if (!Unit.isAtomic(u)) return false;
  return u.model.atomicHierarchy.atoms.label_comp_id.value(u.elements[0]) === 'HOH';
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

async function signature(s: Structure, props: Partial<InteractionsProps>, stats?: InterUnitContactsStats) {
  const ctx = { runtime: SyncRuntimeContext, assetManager: new AssetManager() };
  const { contacts, bridges } = await computeInteractions(ctx, s, props, { stats });
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
    // identity operators
    [
      '7qpd',
      async () => Structure.ofModel(m7qpd),
      'count=94 types=2:6,4:88 ordered=d26d8de5 sorted=26280639 bridges=0',
      'count=200 types=1:10,2:6,4:88,6:70,8:26 ordered=71e8eb91 sorted=ff4f819 bridges=0',
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

  // symmetry 'per-pair-class': searched once per class of symmetry-equivalent unit pairs, line of sight per pair; equal
  // to 'all' up to floating-point ties at the thresholds, which these structures do not have
  it.each(cases)('%s, per pair class', async (_name, structure, defaults, all) => {
    const s = await structure();
    const stats = { pairs: 0, classes: 0, replayed: 0 };
    expect(await signature(s, { symmetry: 'per-pair-class' }, stats)).toBe(defaults);
    expect(await signature(s, { ...AllProviders, symmetry: 'per-pair-class' })).toBe(all);
    expect(stats.pairs).toBeGreaterThan(0);
    expect(stats.classes + stats.replayed).toBeLessThanOrEqual(stats.pairs);
  });

  it('per pair class: replays most pairs of a symmetric cell, computes pairs with covalent inter-unit bonds in full', async () => {
    const cell = await cells(Structure.ofModel(m1cbs), Vec3.create(1, 1, 1));
    const stats = { pairs: 0, classes: 0, replayed: 0 };
    expect(await signature(cell, { symmetry: 'per-pair-class' }, stats)).toBe(await signature(cell, {}));
    expect(stats.replayed).toBeGreaterThan(stats.pairs / 2);
    expect(stats.classes + stats.replayed).toBe(stats.pairs);

    // split chains: units bonded across their boundaries are not replayed
    const split = await cells(splitUnits(Structure.ofModel(m1cbs), 300), Vec3.create(1, 1, 0));
    const splitStats = { pairs: 0, classes: 0, replayed: 0 };
    expect(await signature(split, { ...AllProviders, symmetry: 'per-pair-class' }, splitStats)).toBe(
      await signature(split, AllProviders),
    );
    expect(splitStats.classes + splitStats.replayed).toBeLessThan(splitStats.pairs);
  });

  it('per pair class: contact lookups and loci labels as with all', async () => {
    const cell = await cells(Structure.ofModel(m1cbs), Vec3.create(1, 1, 0));
    const ctx = { runtime: SyncRuntimeContext, assetManager: new AssetManager() };
    const all = await computeInteractions(ctx, cell, AllProviders);
    const perClass = await computeInteractions(ctx, cell, { ...AllProviders, symmetry: 'per-pair-class' });
    const labels = (interactions: Interactions) =>
      interactions.contacts.edges.map((e, i) => {
        const unitA = cell.unitMap.get(e.unitA);
        const unitB = cell.unitMap.get(e.unitB);
        expect(interactions.contacts.getEdgeIndex(e.indexA, e.unitA, e.indexB, e.unitB)).toBe(i);
        return Interactions.Loci(cell, interactions, [{ unitA, indexA: e.indexA, unitB, indexB: e.indexB }]).getLabel();
      });
    expect(labels(perClass)).toEqual(labels(all));
  });

  it('per pair class: a closed C6 assembly gives the same contacts as all', async () => {
    const s = Structure.ofModel(m1cbs);
    // neighbouring copies overlap by their surfaces
    const c6 = cyclic(s, 6, s.boundary.sphere.radius);
    const stats = { pairs: 0, classes: 0, replayed: 0 };
    const all = await signature(c6, AllProviders);
    const index = (id: number) => Math.floor(c6.unitIndexMap.get(id) / s.units.length);
    const { contacts } = await computeInteractions(
      { runtime: SyncRuntimeContext, assetManager: new AssetManager() },
      c6,
      {},
    );
    expect(contacts.edges.filter((e) => index(e.unitA) !== index(e.unitB)).length).toBeGreaterThan(0);
    expect(await signature(c6, { ...AllProviders, symmetry: 'per-pair-class' }, stats)).toBe(all);
    expect(stats.replayed).toBeGreaterThan(0);
  });

  it('per pair class: line of sight is checked per pair (a water blocks one pair of a class)', async () => {
    const s = Structure.ofModel(m1cbs);
    const offset = s.boundary.sphere.radius;
    const c6 = cyclic(s, 6, offset);
    const n = s.units.length;
    const ctx = { runtime: SyncRuntimeContext, assetManager: new AssetManager() };

    // a contact between copy 0 and copy 1, and its feature midpoint
    const { contacts, unitsFeatures } = await computeInteractions(ctx, c6, {});
    const index = (unitId: number) => c6.unitIndexMap.get(unitId);
    const e = contacts.edges.find((e) => index(e.unitA) < n && index(e.unitB) >= n && index(e.unitB) < 2 * n)!;
    expect(e).toBeDefined();
    const [iA, iB] = [index(e.unitA), index(e.unitB)];
    const position = (unitId: number, f: number) => {
      const features = unitsFeatures.get(unitId);
      const p = Vec3.create(features.x[f], features.y[f], features.z[f]);
      return Vec3.transformMat4(p, p, c6.unitMap.get(unitId).conformation.operator.matrix);
    };
    const mid = Vec3.scale(Vec3(), Vec3.add(Vec3(), position(e.unitA, e.indexA), position(e.unitB, e.indexB)), 0.5);

    // one water oxygen moved to that midpoint
    const water = s.units.find(isWater)!;
    const w = water.getChild(SortedArray.ofSingleton(water.elements[0]));
    const wPos = w.conformation.position(w.elements[0], Vec3());
    const blocked = cyclic(s, 6, offset, (builder) =>
      builder.addWithOperator(
        w,
        SymmetryOperator.create('block', Mat4.fromTranslation(Mat4(), Vec3.sub(Vec3(), mid, wPos))),
      ),
    );

    const stats = { pairs: 0, classes: 0, replayed: 0 };
    expect(await signature(blocked, { symmetry: 'per-pair-class' }, stats)).toBe(await signature(blocked, {}));
    expect(stats.replayed).toBeGreaterThan(0);
    const result = await computeInteractions(ctx, blocked, { symmetry: 'per-pair-class' });
    const has = (unitA: number, unitB: number) =>
      result.contacts.getEdgeIndex(e.indexA, blocked.units[unitA].id, e.indexB, blocked.units[unitB].id) >= 0;
    // blocked for copies 0–1, the same contact of copies 1–2 (same class) is still there
    expect(has(iA, iB)).toBe(false);
    expect(has(iA + n, iB + n)).toBe(true);
  });
});
