/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  createRegistries,
  KernelRegistry,
  type LayerTypeEntry,
  MvsRegistry,
  type RegistryKind,
  register,
  SimpleParamsSchema,
  SvsRegistrationError,
  type SvsRegistries,
} from '../../src/format/index.js';
import { registeredNamespaces } from '../../src/format/registry.js';
import { KernelLayerTypes } from '../../src/format/tables/kernel.js';

const KINDS: RegistryKind[] = ['layerTypes', 'themes', 'custom', 'issueCodes'];
const snapshot = (reg: SvsRegistries) => Object.fromEntries(KINDS.map((k) => [k, [...reg.list(k)]]));

const type = (name: string): LayerTypeEntry => ({
  name,
  options: SimpleParamsSchema({}),
  accepts: ['ranges'],
  layouts: ['tracks'],
  axes: ['blocks'],
  sizing: 'row',
  fallback: 'ranges',
  channels: [],
});

describe('registries', () => {
  it("register the kernel and mvs in order, each kind in the interface's order", () => {
    const reg = createRegistries();
    expect(reg.version).toBe(0);
    register(reg, [KernelRegistry, MvsRegistry]);
    expect(reg.version).toBe(1);
    expect(reg.list('layerTypes').map((t) => t.name)).toEqual(KernelLayerTypes.map((t) => t.name));
    expect(reg.list('themes').filter((t) => t.params.color)).toHaveLength(12);
    expect(reg.list('themes').filter((t) => t.params.glyph)).toHaveLength(4);
    expect(reg.get('custom', 'mvs/unit')?.namespace).toBe('mvs');
    expect(reg.list('layerTypes')).toBe(reg.list('layerTypes'));
    expect([...registeredNamespaces(reg)].sort()).toEqual(['mvs', 'svs']);
  });

  it('the same entry object twice is not a conflict', () => {
    const reg = createRegistries();
    register(reg, KernelRegistry);
    const before = snapshot(reg);
    const undo = register(reg, [KernelRegistry, KernelRegistry]);
    expect(snapshot(reg)).toEqual(before);
    undo();
    expect(snapshot(reg)).toEqual(before);
  });

  it('lists every conflict in one error and registers nothing', () => {
    const reg = createRegistries();
    register(reg, KernelRegistry);
    const before = snapshot(reg);
    const version = reg.version;
    let message = '';
    let conflicts: readonly string[] = [];
    try {
      register(reg, [
        { layerTypes: [type('org.example:new'), type('letters')] },
        {
          layerTypes: [type('org.example:twice')],
          themes: [{ name: 'heat', params: { color: SimpleParamsSchema({}) }, fallback: null }],
        },
        {
          layerTypes: [type('org.example:twice')],
          custom: [{ namespace: 'svs', kind: 'unit', schema: SimpleParamsSchema({}) }],
        },
        { issueCodes: [{ code: 'svs.mine', severity: 'error' }] },
      ]);
    } catch (e) {
      expect(e).toBeInstanceOf(SvsRegistrationError);
      message = (e as SvsRegistrationError).message;
      conflicts = (e as SvsRegistrationError).conflicts;
    }
    expect(message).toContain('5 conflict(s), nothing registered');
    expect(conflicts).toHaveLength(5);
    expect(conflicts).toContain('layerTypes: "letters" is already registered');
    expect(message).toContain('layerTypes: "letters" is already registered');
    expect(message).toContain('layerTypes: "org.example:twice" is listed twice');
    expect(message).toContain('themes: "heat" is not namespaced');
    expect(message).toContain('custom: "svs" is reserved');
    expect(message).toContain('issueCodes: "svs.mine" is reserved');
    expect(reg.version).toBe(version);
    expect(snapshot(reg)).toEqual(before);
    expect(reg.get('layerTypes', 'org.example:new')).toBeUndefined();
  });

  it('the viewer owns seqv: only the entries it names may register under it', () => {
    const reg = createRegistries();
    const own = {
      issueCodes: [{ code: 'seqv.patch.dropped', severity: 'warning' as const }],
      layerTypes: [type('seqv:map')],
    };
    const ext = { issueCodes: [{ code: 'seqv.mine', severity: 'error' as const }] };
    expect(() => register(reg, [own, ext], { seqv: new Set([own]) })).toThrow('issueCodes: "seqv.mine" is reserved');
    expect(() => register(reg, own)).toThrow('layerTypes: "seqv" is reserved');
    expect(() =>
      register(reg, { custom: [{ namespace: 'svs', kind: 'x', schema: SimpleParamsSchema({}) }] }, { seqv: new Set() }),
    ).toThrow('custom: "svs" is reserved');
    register(reg, own, { seqv: new Set([own]) });
    expect(reg.get('layerTypes', 'seqv:map')).toBe(own.layerTypes[0]);
    expect(registeredNamespaces(reg).has('seqv')).toBe(true);
  });

  it('undo restores the registries, once', () => {
    const reg = createRegistries();
    register(reg, KernelRegistry);
    const before = snapshot(reg);
    const ext = {
      layerTypes: [type('org.example:fingerprint')],
      issueCodes: [{ code: 'org.example.bad', severity: 'warning' as const }],
    };
    const undo = register(reg, [ext, MvsRegistry]);
    expect(reg.get('layerTypes', 'org.example:fingerprint')).toBe(ext.layerTypes[0]);
    expect(registeredNamespaces(reg).has('org.example')).toBe(true);
    undo();
    const version = reg.version;
    expect(snapshot(reg)).toEqual(before);
    expect(reg.get('custom', 'mvs/unit')).toBeUndefined();
    undo();
    expect(reg.version).toBe(version);
  });

  it('a host replaces a built-in by registering a kernel copy without it', () => {
    const reg = createRegistries();
    const mine = type('letters');
    register(reg, [
      { ...KernelRegistry, layerTypes: KernelLayerTypes.filter((t) => t.name !== 'letters') },
      { layerTypes: [mine] },
    ]);
    expect(reg.get('layerTypes', 'letters')).toBe(mine);
  });

  it('refuses registries it did not make', () => {
    const fake = { version: 0, get: () => undefined, list: () => [] } as unknown as SvsRegistries;
    expect(() => register(fake, KernelRegistry)).toThrow(/createRegistries/);
  });
});
