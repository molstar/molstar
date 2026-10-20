/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Structure } from '@molstar/model/model/structure';
import { BondProvider } from '@molstar/model/model/structure/structure/unit/bonds/provider';
import { cbs, expectOk, headlessViewer, until } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';

// `structure.bonds` and the inspector's Bonds field (bond-providers plan §4.1) on a headless viewer.

const structureCell = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures[0].cell;
const typeParams = (v: PrismViewer) => (structureCell(v).transform.params as any).type;
const providerOf = (v: PrismViewer) => (structureCell(v).obj!.data as Structure).bondProvider.provider.name;
const representationCount = (v: PrismViewer) =>
  v.plugin.managers.structure.hierarchy.current.structures[0].components.reduce(
    (n, c) => n + c.representations.length,
    0,
  );

describe('structure.bonds', () => {
  it('starts with the built-in bonds and switches provider, keeping the components', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expect(providerOf(v)).toBe(BondProvider.DefaultName);
    const representations = representationCount(v);

    expectOk(await v.api.structure.bonds({ provider: 'none' }));
    expect(typeParams(v).params.bonds).toEqual({ name: 'none', params: {} });
    expect(providerOf(v)).toBe('none');
    expect(representationCount(v)).toBe(representations);

    const again = await v.api.structure.bonds({ provider: 'none' });
    expect(again).toMatchObject({ ok: true, unchanged: true });

    expectOk(await v.api.structure.bonds({ provider: 'default' }));
    expect(providerOf(v)).toBe(BondProvider.DefaultName);
  });

  it('refuses a provider the plugin lacks, naming the ones it has', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const r = await v.api.structure.bonds({ provider: 'nonsense' });
    expect(r).toMatchObject({ ok: false, code: 'invalid' });
    if (!r.ok) expect(r.hint).toMatch(/default, none, per-unit/);
  });

  it('survives a change of assembly', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.structure.bonds({ provider: 'none' }));
    expectOk(await v.api.structure.assembly({ id: '1' }));
    expect(typeParams(v).name).toBe('assembly');
    expect(typeParams(v).params.bonds).toEqual({ name: 'none', params: {} });
    expect(providerOf(v)).toBe('none');
  });

  it('undoes', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.structure.bonds({ provider: 'none' }));
    await v.api.history.undo();
    expect(providerOf(v)).toBe(BondProvider.DefaultName);
  });
});

describe("the inspector's Bonds field", () => {
  const form = (v: PrismViewer) => {
    const structure = v.views.scene.roots.value[0];
    return v.views.scene.inspect(structure).form!;
  };

  it('offers the providers and edits through structure.bonds', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const f = form(v);
    const bonds = f.params.bonds as PD.Mapped<any>;
    expect(bonds.isHidden).toBeFalsy();
    expect(bonds.select.options.map((o) => o[0])).toEqual(['default', 'none', 'per-unit']);
    expect(f.values.value.bonds).toEqual({ name: 'default', params: {} });

    f.set('bonds', { name: 'none', params: {} });
    await until(() => providerOf(v) === 'none' && v.api.running.value.length === 0);
  });

  it('is hidden when the built-in bonds are the only choice', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    v.plugin.bondProviders.clear();
    expect((form(v).params.bonds as PD.Mapped<any>).isHidden).toBe(true);
  });
});
