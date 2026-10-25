/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import { mapFile } from './volume-fixtures.js';
import type { PrismViewer } from '../../viewer.js';
import type { LookState } from '../look.js';

// `look` (api.md §2, "Look"): the scene defaults, a node's own values, and overrides kept through a change of the
// defaults.

const Plastic = { metalness: 0, roughness: 0.2, bumpiness: 0 };
const Metallic = { metalness: 1, roughness: 0.6, bumpiness: 0 };

const steps = (v: PrismViewer) => v.managers.history.steps.value.length;
const reprsOf = (v: PrismViewer, label: string) =>
  v.plugin.managers.structure.hierarchy.current.structures
    .flatMap((s) => s.components)
    .filter((c) => (c.cell.params?.values as any).label === label)
    .flatMap((c) => c.representations);
const params = (v: PrismViewer, label: string) => reprsOf(v, label)[0].cell.transform.params as any;
const others = (v: PrismViewer, label: string) =>
  v.plugin.managers.structure.hierarchy.current.structures
    .flatMap((s) => s.components)
    .filter((c) => (c.cell.params?.values as any).label !== label)
    .flatMap((c) => c.representations);

async function withLigand() {
  const v = await headlessViewer({ files: [cbs()] });
  expectOk(await v.api.component.add({ target: 'ligand', label: 'RA', representation: 'ball-and-stick' }));
  return v;
}

describe('look.set', () => {
  it('sets the scene defaults and the representations, one step; the same values are already so', async () => {
    const v = await withLigand();
    const manager = v.plugin.managers.structure.component;
    const before = steps(v);
    expectOk(await v.api.look.set({ ignoreLight: true, quality: 'high' }));
    expect(manager.state.options).toMatchObject({ ignoreLight: true, visualQuality: 'high' });
    expect(params(v, 'RA').type.params).toMatchObject({ ignoreLight: true, quality: 'high' });
    expect(steps(v)).toBe(before + 1);

    expect(await v.api.look.set({ ignoreLight: true })).toMatchObject({ ok: true, unchanged: true });
    expect(await v.api.look.set({})).toMatchObject({ ok: true, unchanged: true });
    expect(steps(v)).toBe(before + 1);
    expect(await v.api.look.set({ hydrogens: 'some' as never })).toMatchObject({ ok: false, code: 'invalid' });
  });

  it('runs under its former id, component.options, with the material by its name', async () => {
    const v = await withLigand();
    expectOk(await v.api.run('component.options', { material: Plastic, hydrogens: 'hide-all' }));
    expect(v.plugin.managers.structure.component.state.options).toMatchObject({
      materialStyle: Plastic,
      hydrogens: 'hide-all',
    });
  });

  it('gives a representation added later the defaults', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.look.set({ ignoreLight: true }));
    expectOk(await v.api.component.add({ target: 'ligand', label: 'RA' }));
    expect(params(v, 'RA').type.params.ignoreLight).toBe(true);
  });

  it('sets the representations under a node, and keeps them through a change of the default', async () => {
    const v = await withLigand();
    const ra = reprsOf(v, 'RA')[0].component.cell.transform.ref;
    const r = await v.api.look.set({ node: ra, material: Metallic });
    expectOk(r);
    expect(r.did).toMatch(/material of Ball & Stick of RA/);
    expect(params(v, 'RA').type.params.material).toEqual(Metallic);

    // the default changes: RA keeps its own, the rest follow, and the result says so
    const d = await v.api.look.set({ material: Plastic });
    expectOk(d);
    expect(d.did).toMatch(/keeps its own/);
    expect(params(v, 'RA').type.params.material).toEqual(Metallic);
    for (const o of others(v, 'RA')) expect((o.cell.transform.params as any).type.params.material).toEqual(Plastic);
  });

  it('refuses the scene-only values on a node', async () => {
    const v = await withLigand();
    const ra = reprsOf(v, 'RA')[0].component.cell.transform.ref;
    const clip = v.plugin.managers.structure.component.state.options.clipObjects;
    expect(await v.api.look.set({ node: ra, clip })).toMatchObject({ ok: false, code: 'invalid' });
  });
});

describe('look.get and look.reset', () => {
  it('reads the overrides and takes them back to the defaults', async () => {
    const v = await withLigand();
    const ra = reprsOf(v, 'RA')[0].component.cell.transform.ref;
    expectOk(await v.api.look.set({ node: ra, material: Metallic, ignoreLight: true }));

    const got = await v.api.look.get({});
    expectOk(got);
    const state = got.value as LookState;
    expect(state.overrides).toHaveLength(1);
    expect(state.overrides[0].values).toEqual({ ignoreLight: true, material: Metallic });
    expect(got.did).toMatch(/has its own ignore light on, material Metallic/);

    // one value back
    expectOk(await v.api.look.reset({ node: ra, keys: ['ignoreLight'] }));
    expect(params(v, 'RA').type.params.ignoreLight).toBe(false);
    expect(params(v, 'RA').type.params.material).toEqual(Metallic);

    // the rest, from the scene (Apply to all)
    expectOk(await v.api.look.reset({ keys: ['material'] }));
    expect(await v.api.look.get({})).toMatchObject({ value: { overrides: [] } });
    expect(await v.api.look.reset({})).toMatchObject({ ok: true, unchanged: true });
  });
});

describe('look by node on other representations', () => {
  it('sets a map isosurface’s light and material: its own, no scene default to follow', async () => {
    const v = await headlessViewer({ files: [mapFile()] });
    const iso = [...v.plugin.state.data.cells.values()].find(
      (c) => c.obj && (c.transform.params as any)?.type?.name === 'isosurface',
    )!;
    expect((iso.transform.params as any).type.params.ignoreLight).toBe(false);
    expectOk(await v.api.look.set({ node: iso.transform.ref, ignoreLight: true }));
    expect((v.plugin.state.data.cells.get(iso.transform.ref)!.transform.params as any).type.params.ignoreLight).toBe(
      true,
    );
    expect(await v.api.look.set({ node: iso.transform.ref, ignoreLight: true })).toMatchObject({ unchanged: true });
    // nothing of a map follows the defaults: nothing to reset, and said so for its node
    expect(await v.api.look.reset({})).toMatchObject({ unchanged: true });
    const r = await v.api.look.reset({ node: iso.transform.ref });
    expect(r).toMatchObject({ ok: false, code: 'invalid' });
    expect(!r.ok && r.message).toMatch(/follows no scene defaults/);
  });
});
