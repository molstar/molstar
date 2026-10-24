/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ShellModel } from '@molstar/prism-ui';
import { describe, expect, it } from 'vitest';
import { sequenceCommands } from '../sequence/commands.js';
import { cbs, crambin, expectOk, headlessViewer, until } from './harness.js';

// The `sequence` area's operations (api/sequence.ts, api.md §2, sequence.md §6), headless.

/** The built document's tracks, as `sequence.get` lists them. */
const Tracks = [
  { ref: 'ruler', name: 'Ruler', hidden: false, removable: true },
  { ref: 'letters', name: 'Residues', hidden: false, removable: true },
  { ref: 'ss-track', name: 'Secondary structure', hidden: false, removable: true },
];

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

  it("sets the tracks' rail as a step: the navigator and its height, already so, and Undo restores it", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const viewer = (await v.views.sequence!.load()).viewers.get();
    expect(viewer.rail.value).toEqual({ navigator: 'ruler' });
    const steps = v.managers.history.steps.value.length;

    const r = await v.api.sequence!.rail({ navigator: 'strip', height: 40 });
    expectOk(r);
    expect(viewer.rail.value).toEqual({ navigator: 'strip', height: 40 });
    expect(r.did).toBe("The tracks' rail is the strip, 40 px.");
    expect(v.managers.history.steps.value).toHaveLength(steps + 1);

    const again = await v.api.sequence!.rail({ navigator: 'strip', height: 40 });
    expect(again).toMatchObject({ ok: true, unchanged: true });
    expect(v.managers.history.steps.value).toHaveLength(steps + 1);

    // A value left out is kept: the height alone, the navigator alone (which keeps its own height).
    expectOk(await v.api.sequence!.rail({ height: 50 }));
    expect(viewer.rail.value).toEqual({ navigator: 'strip', height: 50 });
    expect(await v.api.sequence!.rail({ navigator: 'strip' })).toMatchObject({ ok: true, unchanged: true });
    // Another navigator starts at its own height.
    expectOk(await v.api.sequence!.rail({ navigator: 'band' }));
    expect(viewer.rail.value).toEqual({ navigator: 'band' });
    expect((await v.api.sequence!.get()) as any).toMatchObject({ value: { rail: { navigator: 'band' } } });

    expectOk(await v.api.history.undo());
    expect(viewer.rail.value).toEqual({ navigator: 'strip', height: 50 });
  });

  it("the sequence's rail commands change the rail through the API, so it is a recorded step", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const viewer = (await v.views.sequence!.load()).viewers.get();
    const shell = new ShellModel();
    shell.registerTabs([
      {
        id: 'sequence',
        title: 'Sequence',
        icon: 'sequence',
        dock: 'top',
        frame: 'none',
        render: { mount: () => () => {} },
      } as never,
    ]);
    const commands = sequenceCommands(viewer, shell, v.api);
    const runs: string[] = [];
    v.api.events.finished.subscribe((r) => runs.push(r.operation));
    const steps = v.managers.history.steps.value.length;

    await commands.find((c) => c.id === 'sequence.rail.band')!.run({}, {}, new AbortController().signal);
    await until(() => viewer.rail.value.navigator === 'band' && runs.includes('sequence.rail'));
    expect(runs).toContain('sequence.rail');
    expect(v.managers.history.steps.value).toHaveLength(steps + 1);

    expectOk(await v.api.history.undo());
    expect(viewer.rail.value).toEqual({ navigator: 'ruler' });
  });

  it('refuses a rail call that changes nothing, a height for a rail without a row, and a navigator it does not know', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expect(await v.api.sequence!.rail({})).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.sequence!.rail({ height: 30 })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.sequence!.rail({ navigator: 'none', height: 30 })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.sequence!.rail({ navigator: 'wide' as never })).toMatchObject({ ok: false, code: 'invalid' });
    const viewer = (await v.views.sequence!.load()).viewers.get();
    expect(viewer.rail.value).toEqual({ navigator: 'ruler' });
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
    await until(() => viewer.tracks.value.length === 3);
    expect(await get()).toMatchObject({
      structure: '1CBS',
      followed: 'automatic',
      style: 'tracks',
      colors: 'prism:residues',
      rail: { navigator: 'ruler' },
      tracks: Tracks,
      addable: [],
    });
    // Its chains, as `sequence.chains` takes them, and no filter on them.
    const first = await get();
    expect(first.chains.map((c: { name: string }) => c.name)).toContain('A');
    expect(first).not.toHaveProperty('filter');

    expectOk(await v.api.sequence!.style({ style: 'text' }));
    expectOk(await v.api.sequence!.colors({ scheme: 'taylor' }));
    expectOk(await v.api.sequence!.follow({ structure: '1crn' }));
    await until(() => viewer.structure.value?.title === '1CRN');
    await until(() => viewer.tracks.value.length === 3);
    expect(await get()).toMatchObject({
      structure: '1CRN',
      followed: 'chosen',
      style: 'text',
      colors: 'taylor',
      rail: { navigator: 'ruler' },
      tracks: Tracks,
    });

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
