/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import fs from 'node:fs';
import { describe, expect, it } from 'vitest';
import { parseMol2 } from '@molstar/io/reader/mol2/parser';
import { trajectoryFromMol2 } from '@molstar/model/formats/structure/mol2';
import { Structure, Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import { BondType } from '@molstar/model/model/structure/model/types';
import { BondProviderRegistry } from '@molstar/plugin/state/bonds/registry';
import { BondOrderPerceptionProvider, BondOrders } from '../index.js';
import { bondRows, modelFromPdb, perceived, readFixture, unitOf } from './util.js';

// The provider's wiring (_spec §2): modes, cache, children, registry.

async function adrenalin(): Promise<Model> {
  const name = 'adrenalin.mol2';
  const text = fs.readFileSync(new URL(`../../../../../data/examples/${name}`, import.meta.url), 'utf8');
  const parsed = await parseMol2(text, name).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMol2(parsed.result).run()).representative;
}

const covalentFlags = (u: Unit.Atomic) => {
  const { a, b, edgeProps } = u.bonds;
  const flags: number[] = [];
  for (let i = 0; i < a.length; i++) {
    if (a[i] < b[i] && BondType.isCovalent(edgeProps.flags[i])) flags.push(edgeProps.flags[i]);
  }
  return flags;
};

describe('bond-order perception provider', () => {
  it("auto leaves stated orders (a mol2's bond table) as they are", async () => {
    const model = await adrenalin();
    const builtIn = Structure.ofModel(model).units[0] as Unit.Atomic;
    const auto = perceived(model, 'auto').units[0] as Unit.Atomic;
    expect(bondRows(auto)).toEqual(bondRows(builtIn));
    expect(covalentFlags(auto).some((f) => BondType.is(f, BondType.Flag.Computed))).toBe(false);
  });

  it('forceCompute perceives every covalent bond, stated ones included', async () => {
    const forced = perceived(await adrenalin(), 'forceCompute').units[0] as Unit.Atomic;
    expect(covalentFlags(forced).every((f) => BondType.is(f, BondType.Flag.OrderUnknown))).toBe(true);
    expect(bondRows(forced).filter((r) => r.aromatic).length).toBe(6);
  });

  it('perceived bonds are not carried to another frame; a canonical residue keeps the built-in props', async () => {
    const model = await modelFromPdb(readFixture('5kre-ser124-6wg.pdb'), '5KRE');
    const s = perceived(model);
    expect(unitOf(s, '6WG').bonds.props?.canRemap).toBe(false);
    const builtIn = unitOf(Structure.ofModel(model), 'SER').bonds.props?.canRemap;
    expect(unitOf(s, 'SER').bonds.props?.canRemap).toBe(builtIn);
  });

  it('two structures of one model share the perceived bonds; another mode does not', async () => {
    const model = await modelFromPdb(readFixture('1nco-chr.pdb'), '1NCO');
    const first = unitOf(perceived(model, 'auto'), 'CHR').bonds;
    expect(unitOf(perceived(model, 'auto'), 'CHR').bonds).toBe(first);
    expect(unitOf(perceived(model, 'forceCompute'), 'CHR').bonds).not.toBe(first);
  });

  it("one instance per mode, with the root unit's bonds for a child", () => {
    const auto = BondOrderPerceptionProvider.create({ mode: 'auto' });
    expect(BondOrderPerceptionProvider.create({ mode: 'auto' })).toBe(auto);
    expect(auto.cacheKey).toBe('auto');
    expect(auto.childBonds).toBe('root');
  });

  it('BondOrders registers the provider; a choice resolves to its instance', () => {
    expect(BondOrders.structure?.bondProviders).toEqual([BondOrderPerceptionProvider]);
    const registry = new BondProviderRegistry();
    registry.add(BondOrderPerceptionProvider);
    const instance = registry.resolve({ name: 'bond-order-perception', params: { mode: 'forceCompute' } });
    expect(instance).toBe(BondOrderPerceptionProvider.create({ mode: 'forceCompute' }));
  });
});
