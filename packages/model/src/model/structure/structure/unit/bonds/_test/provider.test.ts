/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { SortedArray } from '@molstar/core/data/int';
import { IntAdjacencyGraph } from '@molstar/core/math/graph';
import { Task } from '@molstar/core/task';
import { UUID } from '@molstar/core/util/uuid';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import { ArrayTrajectory } from '../../../../trajectory.js';
import type { ElementIndex } from '../../../../model/indexing.js';
import { StructureSymmetry } from '../../../symmetry.js';
import type { IntraUnitBonds } from '../data.js';
import { computeIntraUnitBonds } from '../intra-compute.js';
import { BondProvider } from '../provider.js';
import { NoComputedBondsProvider } from '../providers.js';
import { BondType } from '../../../../model/types.js';
import { parseMol2 } from '@molstar/io/reader/mol2/parser';
import { trajectoryFromMol2 } from '@molstar/model/formats/structure/mol2';

// Bond providers in the model layer (bond-providers plan §2): the instance a unit and a structure hold, how it passes
// to children, copies, assemblies and frames, the per-key intra-unit cache, and that the default is today's bonds.

const examples = path.resolve(import.meta.dirname, '../../../../../../../../../data/examples');

async function loadCif(name: string): Promise<Model> {
  const parsed = await CIF.parse(fs.readFileSync(path.join(examples, name), 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
}

async function loadMol2(name: string): Promise<Model> {
  const parsed = await parseMol2(fs.readFileSync(path.join(examples, name), 'utf8'), name).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMol2(parsed.result).run()).representative;
}

/** The bonds of a graph as sorted `a-b:order:flags` rows, each once. */
function intraRows(bonds: IntraUnitBonds) {
  const rows: string[] = [];
  const { a, b, edgeProps } = bonds;
  for (let i = 0; i < a.length; i++) {
    if (a[i] < b[i]) rows.push(`${a[i]}-${b[i]}:${edgeProps.order[i]}:${edgeProps.flags[i]}`);
  }
  return rows.sort();
}

function interRows(s: Structure) {
  return s.interUnitBonds.edges
    .filter((e) => e.unitA < e.unitB)
    .map((e) => `${e.unitA}.${e.indexA}-${e.unitB}.${e.indexB}:${e.props.order}:${e.props.flag}`)
    .sort();
}

/** A model that shares no cached bonds with `model`: its own `_dynamicPropertyData`, with the format's (index pairs). */
function fresh(model: Model): Model {
  const data = Object.assign(Object.create(null), model._dynamicPropertyData);
  delete data.ElementSetIntraBondCache;
  delete data.ElementSetIntraBondCacheByKey;
  return { ...model, id: UUID.create22(), _dynamicPropertyData: data };
}

/** The built-in bonds of a unit, marked cacheable or not. */
function withProps(bonds: IntraUnitBonds, cacheable: boolean): IntraUnitBonds {
  const { offset, a, b, edgeCount, edgeProps } = bonds;
  return IntAdjacencyGraph.create(offset, a, b, edgeCount, edgeProps, { canRemap: true, cacheable });
}

/** A provider recording the units and structures it computes bonds for. */
function counting(options: { cacheKey?: string; cacheable?: boolean; childBonds?: 'self' | 'root'; inter?: boolean }) {
  const intraCalls: Unit.Atomic[] = [];
  const interCalls: Structure[] = [];
  const instance: BondProvider.Instance = {
    provider: undefined as any,
    cacheKey: options.cacheKey,
    childBonds: options.childBonds,
    intra(unit, ctx) {
      intraCalls.push(unit);
      return withProps(ctx.builtIn(), !!options.cacheable);
    },
    inter: options.inter
      ? (structure, ctx) => {
          interCalls.push(structure);
          return ctx.builtIn();
        }
      : undefined,
  };
  return { instance, intraCalls, interCalls };
}

const atomic = (s: Structure) => s.units.filter(Unit.isAtomic) as Unit.Atomic[];

const everySecond = (u: Unit) =>
  SortedArray.ofSortedArray<ElementIndex>(Array.from(u.elements).filter((_, i) => i % 2 === 0));

describe('bond providers', async () => {
  const m1cbs = await loadCif('1cbs_updated.cif');

  describe('default', () => {
    it('structures and units get the default instance, with the bonds of the built-in computation', () => {
      const s = Structure.ofModel(fresh(m1cbs));
      expect(s.bondProvider).toBe(BondProvider.DefaultInstance);
      for (const u of atomic(s)) {
        expect(u.bondProvider).toBe(BondProvider.DefaultInstance);
        expect(IntAdjacencyGraph.areEqual(u.bonds, computeIntraUnitBonds(u))).toBe(true);
      }
    });

    it('Unit.create without a provider gives the default', () => {
      const u = Structure.ofModel(m1cbs).units[0];
      const { id, invariantId, chainGroupId, traits, kind, model, conformation, elements } = u;
      const created = Unit.create(id, invariantId, chainGroupId, traits, kind, model, conformation.operator, elements);
      expect((created as Unit.Atomic).bondProvider).toBe(BondProvider.DefaultInstance);
    });
  });

  describe('intra-unit bonds', () => {
    it('are computed by the provider once per unit, and read as the same object', () => {
      const { instance, intraCalls } = counting({});
      const s = Structure.ofModel(fresh(m1cbs), { bondProvider: instance });
      expect(s.bondProvider).toBe(instance);
      const units = atomic(s);
      for (const u of units) {
        expect(u.bondProvider).toBe(instance);
        expect(u.bonds).toBe(u.bonds);
      }
      expect(intraCalls.length).toBe(units.filter((u) => u.elements.length > 0).length);
      expect(new Set(intraCalls).size).toBe(intraCalls.length);
    });

    it('children, copies and remapped units keep the instance', () => {
      const { instance } = counting({});
      const model = fresh(m1cbs);
      const u = Structure.ofModel(model, { bondProvider: instance }).units[0] as Unit.Atomic;
      const child = u.getChild(everySecond(u)) as Unit.Atomic;
      expect(child.bondProvider).toBe(instance);
      expect((child.getCopy(9, 9, 9) as Unit.Atomic).bondProvider).toBe(instance);
      expect((child.applyOperator(9, u.conformation.operator) as Unit.Atomic).bondProvider).toBe(instance);
      expect((u.remapModel(fresh(model), true) as Unit.Atomic).bondProvider).toBe(instance);
    });

    it('a child restricts the bonds of a parent of the same provider', () => {
      const { instance, intraCalls } = counting({});
      const u = Structure.ofModel(fresh(m1cbs), { bondProvider: instance }).units[0] as Unit.Atomic;
      void u.bonds;
      const child = u.getChild(everySecond(u)) as Unit.Atomic;
      void child.bonds;
      expect(intraCalls).toEqual([u]);
    });

    it("childBonds 'self' computes the child's own bonds; 'root' the root unit's, then restricts them", () => {
      const self = counting({ childBonds: 'self' });
      const u1 = Structure.ofModel(fresh(m1cbs), { bondProvider: self.instance }).units[0] as Unit.Atomic;
      const c1 = u1.getChild(everySecond(u1)) as Unit.Atomic;
      void c1.bonds;
      expect(self.intraCalls).toEqual([c1]);
      expect(u1.props.bonds).toBeUndefined();

      const root = counting({ childBonds: 'root' });
      const u2 = Structure.ofModel(fresh(m1cbs), { bondProvider: root.instance }).units[0] as Unit.Atomic;
      const c2 = u2.getChild(everySecond(u2)) as Unit.Atomic;
      void c2.bonds;
      expect(root.intraCalls).toEqual([u2]);
      expect(u2.props.bonds).toBeDefined();
      // the root's bonds restricted, not a computation for the child alone
      expect(c2.bonds.props).toBe(u2.bonds.props);
    });

    it('the cache is per provider key: a model never serves the bonds of one provider to another', () => {
      const model = fresh(m1cbs);
      const a = counting({ cacheKey: 'a', cacheable: true });
      const b = counting({ cacheKey: 'b', cacheable: true });
      const units = (instance: BondProvider.Instance) => atomic(Structure.ofModel(model, { bondProvider: instance }));

      for (const u of units(a.instance)) void u.bonds;
      const computed = a.intraCalls.length;
      expect(computed).toBeGreaterThan(0);

      // another structure of the same provider key: cache hits
      for (const u of units(a.instance)) void u.bonds;
      expect(a.intraCalls.length).toBe(computed);

      // another provider: its own computation
      for (const u of units(b.instance)) void u.bonds;
      expect(b.intraCalls.length).toBe(computed);

      // and the default: the built-in bonds, not a provider's
      for (const u of atomic(Structure.ofModel(model))) {
        expect(IntAdjacencyGraph.areEqual(u.bonds, computeIntraUnitBonds(u))).toBe(true);
        expect(u.bonds.props?.cacheable).not.toBe(true);
      }
    });

    it('no cache key caches nothing', () => {
      const model = fresh(m1cbs);
      const p = counting({ cacheKey: undefined, cacheable: true });
      for (const u of atomic(Structure.ofModel(model, { bondProvider: p.instance }))) void u.bonds;
      const once = p.intraCalls.length;
      for (const u of atomic(Structure.ofModel(model, { bondProvider: p.instance }))) void u.bonds;
      expect(p.intraCalls.length).toBe(2 * once);
    });
  });

  describe('inter-unit bonds', () => {
    it('are computed by the provider of the structure', () => {
      const p = counting({ inter: true });
      const s = Structure.ofModel(fresh(m1cbs), { bondProvider: p.instance });
      expect(s.interUnitBonds).toBe(s.interUnitBonds);
      expect(p.interCalls).toEqual([s]);
    });

    it('a child of the same provider restricts the parent bonds; of another provider computes its own', () => {
      const p = counting({ inter: true });
      const s = Structure.ofModel(fresh(m1cbs), { bondProvider: p.instance });
      void s.interUnitBonds;

      const child = Structure.create(s.units.slice(0, 2), { parent: s });
      expect(child.bondProvider).toBe(p.instance);
      void child.interUnitBonds;
      expect(p.interCalls).toEqual([s]);

      const q = counting({ inter: true });
      const other = Structure.create(s.units.slice(0, 2), { parent: s, bondProvider: q.instance });
      void other.interUnitBonds;
      expect(q.interCalls).toEqual([other]);
    });
  });

  describe("'none': the format's bonds without the distance-based ones", async () => {
    const none = NoComputedBondsProvider.create({});
    const m7qpd = await loadCif('7qpd.fw2.cif');
    const mAce2 = await loadMol2('ace2-hit.mol2');

    it('keeps bond tables (index pairs) as they are', () => {
      const model = fresh(mAce2);
      const byDefault = atomic(Structure.ofModel(model)).map((u) => intraRows(u.bonds));
      const byNone = atomic(Structure.ofModel(model, { bondProvider: none })).map((u) => intraRows(u.bonds));
      expect(byNone).toEqual(byDefault);
      expect(byNone.flat().length).toBeGreaterThan(0);
    });

    it('keeps struct_conn and dictionary bonds of an mmCIF and drops the computed ones, intra and inter', () => {
      const model = fresh(m7qpd);
      const byDefault = Structure.ofModel(model);
      const byNone = Structure.ofModel(model, { bondProvider: none });

      const defaultIntra = atomic(byDefault).flatMap((u) => intraRows(u.bonds).map((r) => `${u.id}|${r}`));
      const noneIntra = atomic(byNone).flatMap((u) => intraRows(u.bonds).map((r) => `${u.id}|${r}`));
      const computed = (row: string) => BondType.is(+row.split(':')[2] as BondType.Flag, BondType.Flag.Computed);
      expect(noneIntra.length).toBeGreaterThan(0);
      expect(noneIntra.some(computed)).toBe(false);
      expect(noneIntra).toEqual(defaultIntra.filter((r) => !computed(r)));
      expect(defaultIntra.some(computed)).toBe(true);

      const defaultInter = interRows(byDefault);
      const noneInter = interRows(byNone);
      expect(noneInter.length).toBeGreaterThan(0);
      expect(noneInter.some(computed)).toBe(false);
      expect(noneInter).toEqual(defaultInter.filter((r) => !computed(r)));
    });

    it('the default bonds are unchanged by the flag being there', () => {
      const model = fresh(m1cbs);
      for (const u of atomic(Structure.ofModel(model))) {
        expect(intraRows(u.bonds)).toEqual(intraRows(computeIntraUnitBonds(u, { noDistance: false })));
      }
    });
  });

  describe('structures derived from a structure', () => {
    it('assemblies, merged units, remapped models and trajectories keep the instance', async () => {
      const { instance } = counting({});
      const model = fresh(m1cbs);
      const s = Structure.ofModel(model, { bondProvider: instance });

      const assembly = await StructureSymmetry.buildAssembly(s, '1').run();
      expect(assembly.bondProvider).toBe(instance);
      for (const u of atomic(assembly)) expect(u.bondProvider).toBe(instance);

      const merged = Structure.mergeUnitsByOperator(s);
      expect(merged.bondProvider).toBe(instance);
      for (const u of atomic(merged)) expect(u.bondProvider).toBe(instance);

      const remapped = s.remapModel(fresh(model));
      expect(remapped.bondProvider).toBe(instance);
      for (const u of atomic(remapped)) expect(u.bondProvider).toBe(instance);

      const trajectory = new ArrayTrajectory([model, fresh(model)]);
      const frames = await Task.create('frames', (ctx) =>
        Structure.ofTrajectory(trajectory, ctx, { bondProvider: instance }),
      ).run();
      expect(frames.bondProvider).toBe(instance);
      for (const u of atomic(frames)) expect(u.bondProvider).toBe(instance);
    });
  });
});
