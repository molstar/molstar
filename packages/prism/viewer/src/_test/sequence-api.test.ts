/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, crambin, expectOk, headlessViewer, until } from './harness.js';

// The `sequence` area's operations (api/sequence.ts, operations.md §2, sequence.md §6), headless.

describe('the sequence API', () => {
  it('colors the letters as a step: set, already so, and Undo restores', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const viewer = (await v.views.sequence!.load()).viewers.get();
    const before = viewer.theme.value.letters;
    const steps = v.managers.history.steps.value.length;

    expectOk(await v.api.sequence!.colors({ scheme: 'zappo' }));
    expect(viewer.theme.value.letters).toBe('zappo');
    expect(v.managers.history.steps.value).toHaveLength(steps + 1);

    const again = await v.api.sequence!.colors({ scheme: 'zappo' });
    expect(again).toMatchObject({ ok: true, unchanged: true });
    expect(v.managers.history.steps.value).toHaveLength(steps + 1);

    expectOk(await v.api.history.undo());
    expect(viewer.theme.value.letters).toBe(before);
    expect(v.managers.history.steps.value).toHaveLength(steps);
  });

  it("turns the letters' colors off with 'none'", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const viewer = (await v.views.sequence!.load()).viewers.get();
    expectOk(await v.api.sequence!.colors({ scheme: 'none' }));
    expect(viewer.theme.value.letters).toBeNull();
    expect(await v.api.sequence!.colors({ scheme: 'none' })).toMatchObject({ ok: true, unchanged: true });
  });

  it('describes the view: the structure it follows, its style and colors, after each change', async () => {
    const v = await headlessViewer({ files: [cbs(), crambin()] });
    const viewer = (await v.views.sequence!.load()).viewers.get();
    const get = async () => {
      const r = await v.api.sequence!.get();
      expectOk(r);
      return r.value!;
    };
    // The viewer takes the followed structure's document asynchronously.
    await until(() => viewer.structure.value?.title === '1CBS');
    expect(await get()).toEqual({ structure: '1CBS', style: 'tracks', colors: 'prism:residues' });

    expectOk(await v.api.sequence!.style({ style: 'text' }));
    expectOk(await v.api.sequence!.colors({ scheme: 'taylor' }));
    expectOk(await v.api.sequence!.follow({ structure: '1crn' }));
    await until(() => viewer.structure.value?.title === '1CRN');
    expect(await get()).toEqual({ structure: '1CRN', style: 'text', colors: 'taylor' });

    expectOk(await v.api.sequence!.colors({ scheme: 'none' }));
    expect((await get()).colors).toBe('none');
  });

  it('says what it shows, in words, and nothing when nothing is loaded', async () => {
    const v = await headlessViewer();
    const r = await v.api.sequence!.get();
    expectOk(r);
    expect(r.did).toBe('No sequence is shown.');
  });

  it('fails to follow a structure that is not loaded, naming those that are', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const r = await v.api.sequence!.follow({ structure: '9xyz' });
    expect(r).toMatchObject({ ok: false, code: 'not-found' });
    if (!r.ok) expect(r.hint).toMatch(/1CBS/);
    // Following stays as it was.
    expect(v.views.sequence!.source.structure.value?.title).toBe('1CBS');
  });

  it('fails to follow when nothing is loaded', async () => {
    const v = await headlessViewer();
    expect(await v.api.sequence!.follow({ structure: '1cbs' })).toMatchObject({ ok: false, code: 'not-found' });
  });

  it('is ambiguous for a name two loaded structures share, and follows either by state ref', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.run('load.files', { files: [cbs()] }, { copy: true }));
    const roots = v.views.scene.roots.value.filter((n) => n.kind === 'structure');
    expect(roots).toHaveLength(2);
    const r = await v.api.sequence!.follow({ structure: '1cbs' });
    expect(r).toMatchObject({ ok: false, code: 'ambiguous' });
    const dataOf = (ref: string) => v.plugin.state.data.cells.get(ref)?.obj?.data;
    const shown = () => v.views.sequence!.source.structure.value?.structure;
    // Following stays as it was: the first.
    expect(shown()).toBe(dataOf(roots[0].id));
    for (const node of [roots[1], roots[0]]) {
      expectOk(await v.api.sequence!.follow({ structure: node.id }));
      await until(() => shown() === dataOf(node.id));
    }
  });
});
