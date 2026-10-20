/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StateSelection } from '@molstar/core/state';
import type { Model } from '@molstar/model/model/structure';
import { Structure, Unit } from '@molstar/model/model/structure';
import { BondProvider } from '@molstar/model/model/structure/structure/unit/bonds/provider';
import { PluginContext } from '@molstar/plugin/context';
import { PluginState } from '@molstar/plugin/state';
import { Mmcif } from '@molstar/plugin/state/formats/trajectory/mmcif';
import { RootStructureDefinition } from '@molstar/plugin/state/helpers/root-structure';
import {
  ModelFromTrajectory,
  StructureFromModel,
  StructureFromTrajectory,
} from '@molstar/plugin/state/transforms/structure/hierarchy';
import { BondProviderRegistry } from '../registry.js';
import {
  NoComputedBondsProvider,
  PerUnitBondProvider,
} from '@molstar/model/model/structure/structure/unit/bonds/providers';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { SortedArray } from '@molstar/core/data/int';

// The plugin's bond providers (bond-providers plan §3, §4): the registry's rules, the param a transform offers and its
// fallback to the built-in bonds, and the instance a structure built through the state tree gets.

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');

/** A provider recording the props it creates instances for; applicable to models whose entry id is not `skip`. */
function testProvider(name: string, created: any[] = []) {
  const provider: BondProvider<{ cutoff: PD.Numeric }> = {
    name,
    label: `Test ${name}`,
    getParams: () => ({ cutoff: PD.Numeric(1.9) }),
    isApplicable: (model) => model.entryId !== 'skip',
    create: (props) => {
      created.push(props);
      return { provider, cacheKey: `${name}:${props.cutoff}`, intra: (_, ctx) => ctx.builtIn() };
    },
  };
  return { provider, created };
}

const model = (entryId: string) => ({ entryId }) as Model;

describe('BondProviderRegistry', () => {
  it("has 'default' always, first, and never as an entry", () => {
    const r = new BondProviderRegistry();
    expect(r.has('default')).toBe(true);
    expect(r.get('default')).toBe(BondProvider.Default);
    expect(r.list).toEqual([]);
    expect(r.options()).toEqual([['default', 'Default']]);
    expect(r.findConflict({ ...BondProvider.Default })).toMatch(/built-in/);
    expect(() => r.add(BondProvider.Default)).toThrow(/built-in/);
  });

  it('counts the same provider, refuses another under its name, removes at zero', () => {
    const r = new BondProviderRegistry();
    const { provider } = testProvider('a');
    r.add(provider);
    r.add(provider);
    expect(r.list).toEqual([provider]);
    expect(() => r.add(testProvider('a').provider)).toThrow(/already registered/);
    r.remove(provider);
    expect(r.has('a')).toBe(true);
    r.remove(provider);
    expect(r.has('a')).toBe(false);
    expect(r.get('a')).toBeUndefined();
    r.remove('unknown');
  });

  it('offers the providers applicable to a model, all of them without one', () => {
    const r = new BondProviderRegistry();
    r.add(testProvider('a').provider);
    expect(r.options().map((o) => o[0])).toEqual(['default', 'a']);
    expect(r.options(model('1abc')).map((o) => o[0])).toEqual(['default', 'a']);
    expect(r.options(model('skip')).map((o) => o[0])).toEqual(['default']);
  });

  it('resolves no choice, the default, unknown and non-applicable names to the default instance', () => {
    const r = new BondProviderRegistry();
    const { provider, created } = testProvider('a');
    r.add(provider);
    expect(r.resolve(undefined)).toBe(BondProvider.DefaultInstance);
    expect(r.resolve({ name: 'default', params: {} })).toBe(BondProvider.DefaultInstance);
    expect(r.resolve({ name: 'missing', params: {} })).toBe(BondProvider.DefaultInstance);
    expect(r.resolve({ name: 'a', params: {} }, model('skip'))).toBe(BondProvider.DefaultInstance);
    expect(created).toEqual([]);

    // params normalized against the provider's: missing ones get their defaults
    const instance = r.resolve({ name: 'a', params: {} }, model('1abc'));
    expect(instance.provider).toBe(provider);
    expect(created).toEqual([{ cutoff: 1.9 }]);
    expect(r.resolve({ name: 'a', params: { cutoff: 2.5 } }).cacheKey).toBe('a:2.5');
  });

  it('the param falls back to the default by normalization: unknown, non-applicable, removed', () => {
    const r = new BondProviderRegistry();
    const { provider } = testProvider('a');
    r.add(provider);
    const normalize = (value: any, m?: Model) =>
      PD.normalizeParams({ bonds: r.getParam(m) }, { bonds: value }, 'all').bonds;

    expect(normalize({ name: 'a', params: { cutoff: 2 } })).toEqual({ name: 'a', params: { cutoff: 2 } });
    expect(normalize({ name: 'missing', params: {} })).toEqual({ name: 'default', params: {} });
    expect(normalize({ name: 'a', params: { cutoff: 2 } }, model('skip'))).toEqual({ name: 'default', params: {} });
    r.remove(provider);
    expect(normalize({ name: 'a', params: { cutoff: 2 } })).toEqual({ name: 'default', params: {} });
  });

  it('builds params without a plugin and without a model', () => {
    const params = RootStructureDefinition.getParams();
    const model = params.type.map('model') as PD.Group<any>;
    expect((model.params.bonds as PD.Mapped<any>).select.options.map((o) => o[0])).toEqual(['default']);
  });
});

describe('bond providers in the plugin', async () => {
  const cif = fs.readFileSync(path.join(examples, '1cbs_updated.cif'), 'utf8');
  const a = testProvider('a');

  async function plugin(withProvider: boolean) {
    const p = new PluginContext({
      behaviors: [],
      registry: [Mmcif, ...(withProvider ? [{ structure: { bondProviders: [a.provider] } }] : [])],
    });
    await p.init();
    const data = await p.builders.data.rawData({ data: cif });
    const trajectory = await p.builders.structure.parseTrajectory(data, 'mmcif');
    return { plugin: p, trajectory };
  }

  async function structureOf(p: PluginContext, trajectory: any, bonds?: any): Promise<Structure> {
    const s = await p
      .build()
      .to(trajectory)
      .apply(ModelFromTrajectory, { modelIndex: 0 })
      .apply(StructureFromModel, { type: { name: 'model', params: bonds ? { bonds } : {} } })
      .commit();
    return s.data!;
  }

  it('the spec registers its providers', async () => {
    const { plugin: p } = await plugin(true);
    expect(p.bondProviders.get('a')).toBe(a.provider);
  });

  it('without a bonds param, and with the default, a structure gets the built-in bonds', async () => {
    const { plugin: p, trajectory } = await plugin(true);
    for (const s of [
      await structureOf(p, trajectory),
      await structureOf(p, trajectory, { name: 'default', params: {} }),
    ]) {
      expect(s.bondProvider).toBe(BondProvider.DefaultInstance);
      for (const u of s.units) if (Unit.isAtomic(u)) expect(u.bondProvider).toBe(BondProvider.DefaultInstance);
    }
  });

  it('a named provider gives the structure and its units its instance', async () => {
    const { plugin: p, trajectory } = await plugin(true);
    const s = await structureOf(p, trajectory, { name: 'a', params: { cutoff: 2.2 } });
    expect(s.bondProvider.provider).toBe(a.provider);
    expect(s.bondProvider.cacheKey).toBe('a:2.2');
    for (const u of s.units) if (Unit.isAtomic(u)) expect(u.bondProvider).toBe(s.bondProvider);
  });

  it('StructureFromTrajectory takes the same choice', async () => {
    const { plugin: p, trajectory } = await plugin(true);
    const s = await p
      .build()
      .to(trajectory)
      .apply(StructureFromTrajectory, { bonds: { name: 'a', params: { cutoff: 2 } } })
      .commit();
    expect(s.data!.bondProvider.provider).toBe(a.provider);
  });

  it('a snapshot naming a provider this plugin lacks loads with the built-in bonds and warns', async () => {
    const source = await plugin(true);
    await structureOf(source.plugin, source.trajectory, { name: 'a', params: { cutoff: 2.2 } });
    const snapshot = JSON.parse(JSON.stringify(source.plugin.state.data.getSnapshot()));

    const { plugin: target } = await plugin(false);
    const warnings: string[] = [];
    (target.log as any).warn = (m: string) => warnings.push(m);
    PluginState.prototype.reportUnregisteredNames.call({ plugin: target } as any, { id: 'id' as any, data: snapshot });
    expect(warnings).toEqual([expect.stringMatching(/'a' is not registered/)]);

    await target.runTask(target.state.data.setSnapshot(snapshot));
    const structures = target.state.data.select(StateSelection.Generators.ofTransformer(StructureFromModel));
    expect(structures.length).toBe(1);
    expect((structures[0].obj!.data as Structure).bondProvider).toBe(BondProvider.DefaultInstance);
  });

  it('a provider removed at run time leaves existing structures as they are; the next update falls back', async () => {
    const { plugin: p, trajectory } = await plugin(true);
    const s = await structureOf(p, trajectory, { name: 'a', params: { cutoff: 2.2 } });
    const [cell] = p.state.data.select(StateSelection.Generators.ofTransformer(StructureFromModel));
    const instance = s.bondProvider;

    p.bondProviders.remove(a.provider);
    expect((cell.obj!.data as Structure).bondProvider).toBe(instance);

    await p
      .build()
      .to(cell.transform.ref)
      .update(StructureFromModel, (old: any) => ({
        type: { ...old.type, params: { ...old.type.params, dynamicBonds: true } },
      }))
      .commit();
    const [updated] = p.state.data.select(StateSelection.Generators.ofTransformer(StructureFromModel));
    expect((updated.obj!.data as Structure).bondProvider).toBe(BondProvider.DefaultInstance);

    // added again: offered again
    p.bondProviders.add(a.provider);
    expect(p.bondProviders.options().map((o) => o[0])).toContain('a');
  });
});

describe("'per-unit'", async () => {
  const cif = fs.readFileSync(path.join(examples, '1cbs_updated.cif'), 'utf8');
  const parsed = await CIF.parse(cif).run();
  if (parsed.isError) throw new Error(parsed.message);
  const m1cbs = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;

  /** A registry with 'none', 'per-unit' and a provider recording the units it computes bonds for. */
  function setup() {
    const calls: Unit.Atomic[] = [];
    const recording: BondProvider<{}> = {
      name: 'rec',
      label: 'Recording',
      getParams: () => ({}),
      create: () => ({
        provider: recording,
        cacheKey: undefined,
        intra: (unit, ctx) => {
          calls.push(unit);
          return ctx.builtIn();
        },
      }),
    };
    const r = new BondProviderRegistry();
    r.add(NoComputedBondsProvider);
    r.add(PerUnitBondProvider);
    r.add(recording);
    return { r, calls };
  }

  const describeUnit = (u: Unit.Atomic) => {
    const e = u.model.atomicHierarchy.index.getEntityFromChain(u.chainIndex[u.elements[0]]);
    return {
      type: u.model.entities.data.type.value(e),
      entity: u.model.entities.data.id.value(e),
      asym: u.model.atomicHierarchy.chains.label_asym_id.value(u.chainIndex[u.elements[0]]),
    };
  };

  function run(match: any) {
    const { r, calls } = setup();
    const choice = { name: 'per-unit', params: { rules: [{ match, provider: { name: 'rec', params: {} } }] } };
    const instance = r.resolve(choice, m1cbs);
    const s = Structure.ofModel(m1cbs, { bondProvider: instance });
    const units = s.units.filter(Unit.isAtomic) as Unit.Atomic[];
    for (const u of units) void u.bonds;
    return { units, calls, instance };
  }

  it('offers the other providers to its rules, not itself', () => {
    const { r } = setup();
    const params = r.get('per-unit')!.getParams(m1cbs, r) as any;
    const provider = params.rules.element.provider as PD.Mapped<any>;
    expect(provider.select.options.map((o: any) => o[0])).toEqual(['default', 'none', 'rec']);
  });

  it('matches by entity type, chain and entity; other units get the built-in bonds', () => {
    for (const [match, pick] of [
      [{ name: 'entity-type', params: ['non-polymer'] }, (d: any) => d.type === 'non-polymer'],
      [{ name: 'label-asym-id', params: 'A' }, (d: any) => d.asym === 'A'],
      [{ name: 'entity-id', params: '1, 2' }, (d: any) => d.entity === '1' || d.entity === '2'],
    ] as const) {
      const { units, calls } = run(match);
      const expected = units.filter((u) => pick(describeUnit(u)));
      expect(expected.length).toBeGreaterThan(0);
      expect(expected.length).toBeLessThan(units.length);
      expect(calls).toEqual(expected);
    }
  });

  it('a child unit gets the provider its root matched', () => {
    const { units, calls } = run({ name: 'entity-type', params: ['polymer'] });
    const polymer = units.find((u) => describeUnit(u).type === 'polymer')!;
    calls.length = 0;
    const child = polymer.getChild(
      SortedArray.ofSortedArray(Array.from(polymer.elements).slice(0, 10)) as any,
    ) as Unit.Atomic;
    // a fresh child of a unit whose bonds exist restricts them; drop them to see the provider chosen
    const orphan = Unit.create(
      child.id,
      child.invariantId,
      child.chainGroupId,
      child.traits,
      child.kind,
      child.model,
      child.conformation.operator,
      child.elements,
      undefined,
      child.bondProvider,
    ) as Unit.Atomic;
    void orphan.bonds;
    expect(calls).toEqual([orphan]);
  });
});
