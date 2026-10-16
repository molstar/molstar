/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from '../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createPrismViewerHeadless } from '../../headless.js';
import type { PrismViewer } from '../../viewer.js';
import { ensemble } from './trajectory-fixtures.js';

// `trajectory` (operations.md §2) on a headless viewer, with a small NMR-like ensemble of three models.

const modelIndex = (v: PrismViewer) =>
  (v.plugin.managers.structure.hierarchy.current.structures[0].model!.cell.transform.params as any).modelIndex;

describe('trajectory', () => {
  it('describes the structures that have several models', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    const r = await v.api.trajectory.get();
    expectOk(r);
    expect(r.value).toMatchObject({ playing: false, structures: [{ count: 3, index: 0 }] });
    expect(r.did).toMatch(/model 1 of 3/);

    const single = await headlessViewer({ files: [cbs()] });
    const none = await single.api.trajectory.get();
    expectOk(none);
    expect(none.value?.structures).toEqual([]);
  });

  it('runs under its former ids too (operations.md §1 "Words"), listed once', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    expectOk(await v.api.run('structure.model', { index: 1 }));
    expect(modelIndex(v)).toBe(1);
    expect(v.api.describe('structure.model')?.definition.id).toBe('trajectory.frame');
    expect(v.api.definitions.value.has('structure.model')).toBe(false);
    const listed = await v.api.run('structure.get', {});
    expect(listed).toMatchObject({ ok: true });
    expectOk(await v.api.structure.get());
  });

  it('shows a frame by index and by delta, without a history step', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    const steps = v.managers.history.steps.value.length;

    const r = await v.api.trajectory.frame({ index: 2 });
    expectOk(r);
    expect(r.unchanged).toBeUndefined();
    expect(modelIndex(v)).toBe(2);
    const got = await v.api.trajectory.get();
    expect(got.ok && got.value?.structures[0].index).toBe(2);

    // Past the end wraps to the start, before the start to the end.
    expectOk(await v.api.trajectory.frame({ by: 1 }));
    expect(modelIndex(v)).toBe(0);
    expectOk(await v.api.trajectory.frame({ by: -1 }));
    expect(modelIndex(v)).toBe(2);

    expect(v.managers.history.steps.value.length).toBe(steps);
  });

  it('structure.model shows a model as viewing too: no history step', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    const steps = v.managers.history.steps.value.length;
    expectOk(await v.api.trajectory.frame({ index: 2 }));
    expect(modelIndex(v)).toBe(2);
    expect(v.managers.history.steps.value.length).toBe(steps);
  });

  it('keeps the earlier steps undoable while it moves', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    expectOk(await v.api.structure.visibility({ visible: false }));
    expectOk(await v.api.trajectory.frame({ index: 1 }));
    // Moving through frames did not make the step before it unsafe to undo.
    const undone = await v.api.history.undo();
    expectOk(undone);
    // And Undo leaves the frame shown (frames are viewing, sequence.md §7.3).
    expect(modelIndex(v)).toBe(1);
  });

  it('says when it is already so, and fails for a frame out of range', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    const same = await v.api.trajectory.frame({ index: 0 });
    expect(same).toMatchObject({ ok: true, unchanged: true });
    const past = await v.api.trajectory.frame({ index: 3 });
    expect(past).toMatchObject({ ok: false, code: 'invalid' });
    if (!past.ok) expect(past.hint).toMatch(/0 to 2/);
    expect(await v.api.trajectory.frame({})).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.trajectory.frame({ index: 1, by: 1 })).toMatchObject({ ok: false, code: 'invalid' });
    expect(modelIndex(v)).toBe(0);
  });

  it('is blocked, with the reason, when no structure has frames; blocked for a named one without', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    // `when` says it (the command menu greys the commands out with the same words).
    expect(v.api.describe('trajectory.frame')?.available).toBe('No structure has several models.');
    const none = await v.api.trajectory.frame({ index: 1 });
    expect(none).toMatchObject({ ok: false, code: 'blocked', message: 'No structure has several models.' });
    const empty = await headlessViewer();
    expect(await empty.api.trajectory.frame({ index: 1 })).toMatchObject({ ok: false, code: 'blocked' });

    const both = await headlessViewer({ files: [cbs(), ensemble()] });
    expect(both.api.describe('trajectory.frame')?.available).toBe(true);
    expect(await both.api.trajectory.frame({ index: 1, structure: '1cbs' })).toMatchObject({
      ok: false,
      code: 'blocked',
    });
  });

  it('asks which structure when several have frames', async () => {
    const v = await headlessViewer({ files: [ensemble(), ensemble(4)] });
    expect(await v.api.trajectory.frame({ index: 1 })).toMatchObject({ ok: false, code: 'ambiguous' });
    const list = (await v.api.trajectory.get()) as any;
    const second = list.value.structures[1];
    expectOk(await v.api.trajectory.frame({ index: 3, structure: second.ref }));
    expect(((await v.api.trajectory.get()) as any).value.structures[1].index).toBe(3);
  });

  it('plays through the frames with the animation, and stops, leaving the history alone', async () => {
    const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
    try {
      expectOk(await v.api.load.files({ files: [ensemble()] }));
      expectOk(await v.api.structure.visibility({ visible: false }));
      const steps = v.managers.history.steps.value.length;

      const r = await v.api.trajectory.play({ mode: 'loop', seconds: 3 });
      expectOk(r);
      expect(((await v.api.trajectory.get()) as any).value.playing).toBe(true);
      expect(await v.api.trajectory.play()).toMatchObject({ ok: true, unchanged: true });

      // There is no render loop here: tick the animation by hand, a second into three.
      const animation = v.plugin.managers.animation;
      await animation.tick(0, true);
      await animation.tick(1500, true);
      await until(() => modelIndex(v) !== 0);
      expect(modelIndex(v)).toBe(1);

      expectOk(await v.api.trajectory.stop());
      expect(await v.api.trajectory.stop()).toMatchObject({ ok: true, unchanged: true });
      expect(((await v.api.trajectory.get()) as any).value.playing).toBe(false);
      expect(v.managers.history.steps.value.length).toBe(steps);
      // The step before is still undoable: the animation was viewing.
      expect(v.managers.history.canUndo.value).toBe(true);
    } finally {
      await v.dispose();
    }
  });

  it('cannot play without a canvas', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    expect(await v.api.trajectory.play()).toMatchObject({ ok: false, code: 'unavailable' });
  });

  it('has a command for each operation', async () => {
    const v = await headlessViewer();
    for (const id of ['trajectory.frame', 'trajectory.play', 'trajectory.stop']) {
      expect(v.managers.commands.get(id)).toBeDefined();
    }
  });
});
