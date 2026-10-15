/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';

// The focus history (`focus.history`, `focus.again`; parity.md §3) on a headless viewer with 1CBS.

async function targets(v: PrismViewer) {
  const r = await v.api.focus.history();
  expectOk(r);
  return r.value!.map((e) => e.target);
}

describe('focus history', () => {
  it('lists the recent foci newest first, as targets focus.set takes', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const empty = await v.api.focus.history();
    expectOk(empty);
    expect(empty.value).toEqual([]);
    expect(v.api.describe('focus.again')?.available).toBe('Nothing was focused yet.');

    expectOk(await v.api.focus.set({ target: 'ligand' }));
    expectOk(await v.api.focus.set({ target: 'A:10-12' }));
    const list = await v.api.focus.history();
    expectOk(list);
    expect(list.value?.map((e) => e.target)).toEqual(['A:10-12', 'A:200']);
    expect(list.value?.[1].selected.residues).toBe(1);
    expect(list.did).toMatch(/^2 recent foci, newest first: 1\. /);

    // A listed target focuses the same again.
    const ligand = list.value![1];
    const back = await v.api.focus.set({ target: ligand.target, structure: ligand.structure });
    expectOk(back);
    expect(await targets(v)).toEqual(['A:200', 'A:10-12']);
  });

  it('focuses one again by its place, with choices in its form and no history step', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.focus.set({ target: 'ligand' }));
    expectOk(await v.api.focus.set({ target: 'A:10-12' }));
    const steps = v.api.history.steps.value.length;

    const params = v.api.describe('focus.again')?.params;
    expect(params?.entry.type).toBe('select');
    expect(params?.entry.defaultValue).toBe('2');
    expect((params?.entry as unknown as { options: [string, string][] }).options.map((o) => o[0])).toEqual(['1', '2']);

    const r = await v.api.focus.again({ entry: '2' });
    expectOk(r);
    expect(r.changed).toEqual({ view: true });
    const focused = await v.api.focus.get();
    expectOk(focused);
    expect(focused.value?.selected.residues).toBe(1);
    expect(await targets(v)).toEqual(['A:200', 'A:10-12']);
    expect(v.api.history.steps.value.length).toBe(steps);

    const same = await v.api.focus.again({ entry: '1' });
    expectOk(same);
    expect(same.unchanged).toBe(true);
    expect(await v.api.focus.again({ entry: '7' })).toMatchObject({ ok: false });
    expect(v.managers.commands.get('focus.again')).toBeDefined();
  });

  it('forgets the foci of a removed structure', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.focus.set({ target: 'ligand' }));
    expectOk(await v.api.session.clear());
    const list = await v.api.focus.history();
    expectOk(list);
    expect(list.value).toEqual([]);
  });
});
