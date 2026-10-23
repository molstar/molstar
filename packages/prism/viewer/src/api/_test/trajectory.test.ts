/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from '../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createPrismViewerHeadless } from '../../headless.js';
import type { PrismViewer } from '../../viewer.js';
import { ensemble } from './trajectory-fixtures.js';
import { isPlayingFrames, playingStructure } from '../../sources/trajectory.js';
import { createTrajectoryControls } from '../../ui/trajectory.js';

// `trajectory` (api.md §2) on a headless viewer, with a small NMR-like ensemble of three models.

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

  it('writes no record entry for a frame shown (a scrub seeks once a frame), and no history step', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    await v.record.flush();
    const before = (await v.record.read(50)).items.length;
    const steps = v.managers.history.steps.value.length;
    for (const by of [1, 1, -1]) expectOk(await v.api.trajectory.frame({ by }));
    await v.record.flush();
    expect((await v.record.read(50)).items.length).toBe(before);
    expect(v.managers.history.steps.value.length).toBe(steps);
    expect(modelIndex(v)).toBe(1);
  });

  it('logs no transform timing line for a frame shown (perf-audit-2 10/P14)', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    const lines: string[] = [];
    const sub = v.plugin.events.log.subscribe((e) => lines.push(String(e.message)));
    for (const by of [1, 1, -1]) expectOk(await v.api.trajectory.frame({ by }));
    sub.unsubscribe();
    expect(modelIndex(v)).toBe(1);
    expect(lines.filter((l) => /^Updated /.test(l))).toEqual([]);
  });

  it('runs under its former ids too (api.md §1 "Words"), listed once', async () => {
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

  it('plays through the frames with the animation, and stops, leaving the history alone: no checkpoint', async () => {
    const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
    try {
      expectOk(await v.api.load.files({ files: [ensemble()] }));
      expectOk(await v.api.structure.visibility({ visible: false }));
      const steps = v.managers.history.steps.value.length;
      // A checkpoint is the plugin's snapshot (sources/checkpoint.ts): none while it plays (sequence.md §7.3).
      const checkpoints = vi.spyOn(v.plugin.state, 'getSnapshot');

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
      expect(checkpoints).not.toHaveBeenCalled();
      // The step before is still undoable: the animation was viewing. Undo takes back that step and keeps the frame.
      expect(v.managers.history.canUndo.value).toBe(true);
      expectOk(await v.api.history.undo());
      expect(modelIndex(v)).toBe(1);
      expect(v.plugin.managers.structure.hierarchy.current.structures[0].cell.state.isHidden).toBeFalsy();
    } finally {
      await v.dispose();
    }
  });

  it("plays one structure's frames when named, leaving another trajectory where it is", async () => {
    const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
    try {
      expectOk(await v.api.load.files({ files: [ensemble(), ensemble(4)] }));
      const [a, b] = v.plugin.managers.structure.hierarchy.current.structures;
      const index = (s: typeof a) =>
        (v.plugin.state.data.cells.get(s.model!.cell.transform.ref)!.transform.params as { modelIndex: number })
          .modelIndex;
      const r = await v.api.trajectory.play({ structure: b.cell.transform.ref, mode: 'loop', seconds: 3 });
      expectOk(r);
      expect(JSON.stringify(r)).toMatch(/models/);
      expect(isPlayingFrames(v.plugin, b.model!.cell.transform.ref)).toBe(true);
      expect(isPlayingFrames(v.plugin, a.model!.cell.transform.ref)).toBe(false);
      expect(playingStructure(v.plugin)).toBe(b.cell.transform.ref);
      const animation = v.plugin.managers.animation;
      await animation.tick(0, true);
      await animation.tick(1500, true);
      await until(() => index(b) !== 0);
      expect(index(a)).toBe(0);
      // Playing the other stops this one (one plays at a time).
      expectOk(await v.api.trajectory.play({ structure: a.cell.transform.ref }));
      expect(playingStructure(v.plugin)).toBe(a.cell.transform.ref);
      expect(isPlayingFrames(v.plugin, b.model!.cell.transform.ref)).toBe(false);
      expectOk(await v.api.trajectory.stop());
    } finally {
      await v.dispose();
    }
  });

  it('marks every row playing on a play of all, and plays all while one plays (review fixes)', async () => {
    const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
    try {
      expectOk(await v.api.load.files({ files: [ensemble(), ensemble(4)] }));
      const [a, b] = v.views.scene.roots.value.map((n) => n.trajectory!);
      const off = [a.playing.listen?.(() => {}), b.playing.listen?.(() => {})];
      expectOk(await v.api.trajectory.play({ structure: v.views.scene.roots.value[0].id }));
      await until(() => a.playing.value && !b.playing.value);
      // All while one plays: all of them (not "already").
      const r = await v.api.trajectory.play();
      expect(r).toMatchObject({ ok: true });
      expect((r as { unchanged?: boolean }).unchanged).toBeFalsy();
      await until(() => a.playing.value && b.playing.value);
      // A row's Pause stops it (every one: they play together).
      b.togglePlay();
      await until(() => !a.playing.value && !b.playing.value);
      for (const o of off) o?.();
    } finally {
      await v.dispose();
    }
  });

  it('adds up quick steps from a row, from the frame last asked for (review fix)', async () => {
    const v = await headlessViewer({ files: [ensemble(6)] });
    const t = v.views.scene.roots.value[0].trajectory!;
    const off = t.index.listen?.(() => {});
    t.step(1);
    t.step(1);
    t.step(1);
    await until(() => t.index.value === 3, 10000);
    await new Promise((r) => setTimeout(r, 200));
    expect(t.index.value).toBe(3);
    off?.();
  });

  it('keeps the canvas bar on the structure that played last, and knows when its row is on screen', async () => {
    const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
    try {
      expectOk(await v.api.load.files({ files: [ensemble(), ensemble(4)] }));
      const [ra, rb] = v.views.scene.roots.value;
      const bar = createTrajectoryControls(v.plugin, v.api);
      const off = [bar.model.label.listen?.(() => {}), bar.rowShown.listen?.(() => {})];
      expectOk(await v.api.trajectory.play({ structure: rb.id }));
      await until(() => bar.model.count.value === 4);
      expectOk(await v.api.trajectory.stop());
      await new Promise((r) => setTimeout(r, 50));
      // Paused: still the one that played (Play plays it again), not the first.
      expect(bar.model.count.value).toBe(4);
      expect(bar.rowShown.value).toBe(false);
      rb.trajectory!.markShown(true);
      rb.trajectory!.markShown(true);
      expect(bar.rowShown.value).toBe(true);
      rb.trajectory!.markShown(false);
      expect(bar.rowShown.value).toBe(true);
      rb.trajectory!.markShown(false);
      expect(bar.rowShown.value).toBe(false);
      // Another structure's row on screen is not this one's.
      ra.trajectory!.markShown(true);
      expect(bar.rowShown.value).toBe(false);
      ra.trajectory!.markShown(false);
      for (const o of off) o?.();
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
