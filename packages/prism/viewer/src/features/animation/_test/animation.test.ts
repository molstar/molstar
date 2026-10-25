/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { ensemble } from '../../../api/_test/trajectory-fixtures.js';
import { cbs, expectOk, headlessViewer } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import type { PrismViewer } from '../../../viewer.js';
import { AnimationFeature } from '../index.js';
import { wiggleCells, wiggleOf } from '../source.js';

// The `animation` feature (api.md §2): the plugin's animations and the wiggle, on a headless viewer. The plugin
// ticks its animations on the canvas's frames, so a test ticks by hand; the camera ones read a stand-in canvas.

// The plugin's ids, and the names `animation.play` takes and `animation.get` gives back (the id without its namespace).
const names = {
  spin: 'built-in.animate-camera-spin',
  models: 'built-in.animate-model-index',
  snapshots: 'built-in.animate-state-snapshots',
};
const short = {
  spin: 'animate-camera-spin',
  models: 'animate-model-index',
  snapshots: 'animate-state-snapshots',
};

/** A viewer whose host says it has a canvas. */
async function withCanvas(files: File[] = []): Promise<PrismViewer> {
  const v = await createPrismViewerHeadless({
    spec: {
      ...DefaultPrismViewerSpec(),
      extensions: [...(DefaultPrismViewerSpec().extensions ?? []), AnimationFeature].filter(
        (e, i, all) => all.findIndex((x) => x.id === e.id) === i,
      ),
    },
    modules: ['canvas', 'files'],
  });
  for (const f of files) expectOk(await v.api.load.files({ files: [f] }));
  return v;
}

const modelIndex = (v: PrismViewer) =>
  (v.plugin.managers.structure.hierarchy.current.structures[0].model!.cell.transform.params as any).modelIndex;

describe('animation', () => {
  it("lists the animations the plugin registers, which are the operation's choice too", async () => {
    const v = await headlessViewer({ extensions: [AnimationFeature] });
    for (const id of ['animation.get', 'animation.play', 'animation.stop', 'animation.wiggle']) {
      expect(v.api.describe(id)).toBeDefined();
    }
    const r = await v.api.animation!.get();
    expectOk(r);
    // The names are the short ones people pass back to `animation.play`; the ids are the plugin's.
    const listed = r.value!.animations.map((a) => a.name);
    expect(listed).toEqual(expect.arrayContaining([short.spin, short.models, short.snapshots]));
    expect(r.value!.animations.map((a) => a.id)).toEqual(expect.arrayContaining([names.spin, names.models]));
    // Nothing to animate yet: the trajectory's says why.
    expect(r.value!.animations.find((a) => a.name === short.models)?.available).not.toBe(true);
    expect(r.value!.playing).toBeUndefined();
    expect(r.did).toContain(`${short.spin}: `);
    expect(r.did).not.toContain('built-in.');
    // The play operation's name says which are registered.
    const play = v.api.describe('animation.play')!;
    expect(play.params.name.isOptional).toBe(true);
    expect(play.params.name.description).toContain(`${short.models} (Animate Trajectory)`);
  });

  it('says durations as "Duration (ms)", and listens to the plugin only while watched', async () => {
    const v = await headlessViewer({ files: [cbs()], extensions: [AnimationFeature] });
    const spin = v.api.animation!.paramsOf(names.spin)!;
    expect(spin.durationInMs.label).toBe('Duration (ms)');
    const changed = v.plugin.state.data.events.changed as unknown as { observers: unknown[] };
    const idle = changed.observers.length;
    // Reading the state (the commands do) is a read, not a subscription.
    expect(v.api.animation!.state.value.animations.length).toBeGreaterThan(0);
    expect(changed.observers.length).toBe(idle);
    const watching = v.api.animation!.state.subscribe();
    expect(changed.observers.length).toBeGreaterThan(idle);
    watching.unsubscribe();
    expect(changed.observers.length).toBe(idle);
  });

  it('has the commands: play, stop, wiggle, and the card', async () => {
    const v = await headlessViewer({ extensions: [AnimationFeature] });
    for (const id of [
      'animation.play',
      'animation.stop',
      'animation.wiggle',
      'animation.wiggle.off',
      'animation.open',
    ]) {
      expect(v.managers.commands.get(id), id).toBeDefined();
    }
  });

  it('has one play command whose form lists the animations registered when it is shown', async () => {
    const v = await headlessViewer({ extensions: [AnimationFeature] });
    const play = v.managers.commands.get('animation.play')!;
    const options = () =>
      ((play.params!({ empty: false }, {}) as any).name.options as [string, string][]).map((o) => o[0]);
    expect(options()).toEqual(expect.arrayContaining([short.spin, short.models]));
    // One registered later is in the form, with no command made for it.
    const later = {
      ...v.plugin.managers.animation.animations[0],
      name: 'ext.animate-later',
      display: { name: 'Later' },
    };
    v.plugin.managers.animation.register(later as never);
    expect(options()).toContain('animate-later');
    expect(play.defaults!({ empty: false })).toMatchObject({ name: expect.any(String) });
  });

  it('cannot play without a canvas', async () => {
    const v = await headlessViewer({ extensions: [AnimationFeature] });
    expect(await v.api.animation!.play({ name: names.spin })).toMatchObject({ ok: false, code: 'unavailable' });
  });

  it('refuses an animation that is not there, values it does not take, and one with nothing to animate', async () => {
    const v = await withCanvas([cbs()]);
    try {
      expect(await v.api.animation!.play({ name: 'nope' })).toMatchObject({ ok: false, code: 'not-found' });
      expect(await v.api.animation!.play({ name: names.spin, params: { nonsense: 1 } })).toMatchObject({
        ok: false,
        code: 'invalid',
      });
      expect(await v.api.animation!.play({ name: names.spin, params: { speed: 99 } })).toMatchObject({
        ok: false,
        code: 'invalid',
      });
      const blocked = await v.api.animation!.play({ name: names.models });
      expect(blocked).toMatchObject({ ok: false, code: 'blocked' });
      if (!blocked.ok) expect(blocked.message).toMatch(/No trajectory/);
    } finally {
      await v.dispose();
    }
  });

  it('plays the camera spin on a canvas with its params, says when it is already so, and stops', async () => {
    const v = await withCanvas([cbs()]);
    try {
      const reset = vi.fn();
      const snapshot = { position: [0, 0, 10], target: [0, 0, 0], up: [0, 1, 0], radiusMax: 10 };
      Object.defineProperty(v.plugin, 'canvas3d', {
        value: { camera: { getSnapshot: () => snapshot }, requestCameraReset: reset, dispose() {} },
        configurable: true,
      });
      const steps = v.managers.history.steps.value.length;

      const r = await v.api.animation!.play({ name: names.spin, params: { durationInMs: 1000, speed: 2 } });
      expectOk(r);
      expect(v.plugin.managers.animation.current?.paramValues).toMatchObject({ durationInMs: 1000, speed: 2 });
      expect(v.api.animation!.state.value.playing).toBe(short.spin);
      expect(await v.api.animation!.play({ name: names.spin, params: { durationInMs: 1000, speed: 2 } })).toMatchObject(
        {
          ok: true,
          unchanged: true,
        },
      );

      await v.plugin.managers.animation.tick(0, true);
      await v.plugin.managers.animation.tick(250, true);
      expect(reset).toHaveBeenCalled();

      expectOk(await v.api.animation!.stop());
      expect(await v.api.animation!.stop()).toMatchObject({ ok: true, unchanged: true });
      expect(v.api.animation!.state.value.playing).toBeUndefined();
      expect(v.managers.history.steps.value.length).toBe(steps);
    } finally {
      await v.dispose();
    }
  });

  it('plays the models of an ensemble, which the history does not count', async () => {
    const v = await withCanvas([ensemble()]);
    try {
      expectOk(await v.api.structure.visibility({ visible: false }));
      const steps = v.managers.history.steps.value.length;
      // Params given as the mapped values the animation takes; the rest are its defaults.
      expectOk(
        await v.api.animation!.play({
          name: names.models,
          params: { mode: { name: 'once', params: {} }, duration: { name: 'fixed', params: { durationInS: 3 } } },
        }),
      );
      expect(v.plugin.managers.animation.current?.paramValues.mode).toEqual({
        name: 'once',
        params: { direction: 'forward' },
      });
      await v.plugin.managers.animation.tick(0, true);
      await v.plugin.managers.animation.tick(1500, true);
      expect(modelIndex(v)).toBe(1);
      // The same animation stops as the trajectory's.
      expect((await v.api.trajectory.get()) as any).toMatchObject({ value: { playing: true } });
      expectOk(await v.api.animation!.stop());
      expect(v.managers.history.steps.value.length).toBe(steps);
      expect(v.managers.history.canUndo.value).toBe(true);
    } finally {
      await v.dispose();
    }
  });

  it('does not ask the animations whether they can play on each frame of one playing, while it is watched', async () => {
    const v = await withCanvas([ensemble()]);
    try {
      expectOk(await v.api.structure.visibility({ visible: false }));
      const models = v.plugin.managers.animation.animations.find((a) => a.name === names.models)!;
      const canApply = vi.spyOn(models, 'canApply');
      const watching = v.api.animation!.state.subscribe();
      const frames = () => new Promise((r) => setTimeout(r, 40));
      expectOk(
        await v.api.animation!.play({
          name: names.models,
          params: { mode: { name: 'once', params: {} }, duration: { name: 'fixed', params: { durationInS: 3 } } },
        }),
      );
      await frames();
      expect(v.api.animation!.state.value.playing).toBe(short.models);
      canApply.mockClear();
      for (let t = 0; t <= 1500; t += 250) await v.plugin.managers.animation.tick(t, true);
      await frames();
      expect(modelIndex(v)).toBeGreaterThan(0);
      expect(canApply).not.toHaveBeenCalled();
      // Stopped, it reads again.
      expectOk(await v.api.animation!.stop());
      await frames();
      expect(canApply).toHaveBeenCalled();
      expect(v.api.animation!.state.value.playing).toBeUndefined();
      watching.unsubscribe();
    } finally {
      await v.dispose();
    }
  });
});

describe('animation names', () => {
  /** A canvas stand-in for the camera animations (the plugin reads its camera). */
  const standIn = (v: PrismViewer) =>
    Object.defineProperty(v.plugin, 'canvas3d', {
      value: {
        camera: { getSnapshot: () => ({ position: [0, 0, 10], target: [0, 0, 0], up: [0, 1, 0], radiusMax: 10 }) },
        requestCameraReset: vi.fn(),
        dispose() {},
      },
      configurable: true,
    });
  const registered = (v: PrismViewer, name: string, display: string) => {
    const base = v.plugin.managers.animation.animations.find((a) => a.name === names.models)!;
    v.plugin.managers.animation.register({ ...base, name, display: { name: display } } as never);
  };

  it('takes the id without its namespace, the title in any case, or the full id', async () => {
    const v = await withCanvas([cbs()]);
    standIn(v);
    try {
      for (const name of [short.spin, 'ANIMATE-CAMERA-SPIN', 'Camera Spin', 'camera spin', names.spin]) {
        expectOk(await v.api.animation!.play({ name }));
        expect(v.api.animation!.state.value.playing, name).toBe(short.spin);
        expectOk(await v.api.animation!.stop());
      }
      // One with nothing to animate is found and then blocked, by any of its names.
      for (const name of [short.models, 'Animate Trajectory', names.models]) {
        const r = await v.api.animation!.play({ name });
        expect(r, name).toMatchObject({ ok: false, code: 'blocked' });
        if (!r.ok) expect(r.message).toMatch(/^Animate Trajectory cannot play now/);
      }
    } finally {
      await v.dispose();
    }
  });

  it('says that an unknown name is not found, listing the short names to pass back', async () => {
    const v = await withCanvas();
    try {
      const r = await v.api.animation!.play({ name: 'spin' });
      expect(r).toMatchObject({ ok: false, code: 'not-found' });
      if (!r.ok) {
        expect(r.message).toContain("'spin'");
        expect(r.hint).toContain(`${short.spin} (Camera Spin)`);
        expect(r.hint).not.toContain('built-in.');
      }
    } finally {
      await v.dispose();
    }
  });

  it('is ambiguous for a short name two namespaces share, and the full id tells them apart', async () => {
    const v = await withCanvas();
    try {
      registered(v, 'a.twin', 'Twin A');
      registered(v, 'b.twin', 'Twin B');
      const both = await v.api.animation!.play({ name: 'twin' });
      expect(both).toMatchObject({ ok: false, code: 'ambiguous' });
      if (!both.ok) expect(both.hint).toContain('a.twin, b.twin');
      // A title two animations share is ambiguous too.
      registered(v, 'c.x', 'Same');
      registered(v, 'd.y', 'same');
      expect(await v.api.animation!.play({ name: 'SAME' })).toMatchObject({ ok: false, code: 'ambiguous' });
      expect(await v.api.animation!.play({ name: 'a.twin' })).toMatchObject({ ok: false, code: 'blocked' });
    } finally {
      await v.dispose();
    }
  });

  it('asks which animation when none is given and none plays, listing them, never the first', async () => {
    const v = await withCanvas();
    try {
      const r = await v.api.animation!.play({});
      expect(r).toMatchObject({ ok: false, code: 'invalid' });
      if (!r.ok) {
        expect(r.message).toMatch(/^Which animation\?/);
        expect(r.hint).toContain(short.models);
      }
      expect(v.api.animation!.state.value.playing).toBeUndefined();
    } finally {
      await v.dispose();
    }
  });

  it('plays the one playing when none is given, which is already so', async () => {
    const v = await withCanvas([cbs()]);
    standIn(v);
    try {
      expectOk(await v.api.animation!.play({ name: short.spin }));
      expect(await v.api.animation!.play({})).toMatchObject({ ok: true, unchanged: true });
      // With params, the one playing takes them.
      expectOk(await v.api.animation!.play({ params: { speed: 3 } }));
      expect(v.plugin.managers.animation.current?.paramValues).toMatchObject({ speed: 3 });
      expectOk(await v.api.animation!.stop());
    } finally {
      await v.dispose();
    }
  });
});

describe('wiggle', () => {
  it('wiggles by dynamics, says when it is already so, and undoes', async () => {
    const v = await headlessViewer({ files: [cbs()], extensions: [AnimationFeature] });
    expect(wiggleOf(v.plugin).mode).toBe('off');

    const r = await v.api.animation!.wiggle({ mode: 'dynamics', amplitude: 0.5 });
    expectOk(r);
    expect(wiggleOf(v.plugin)).toMatchObject({ mode: 'dynamics', amplitude: 0.5 });
    // It is on the representations, not just the options.
    const repr = v.plugin.managers.structure.hierarchy.current.structures[0].components[0].representations[0];
    expect((repr.cell.transform.params as any).type.params.animation.wiggleAmplitude).toBe(0.5);
    expect(v.api.animation!.state.value.wiggle.mode).toBe('dynamics');

    expect(await v.api.animation!.wiggle({ mode: 'dynamics', amplitude: 0.5 })).toMatchObject({
      ok: true,
      unchanged: true,
    });
    const steps = v.managers.history.steps.value.length;
    expectOk(await v.api.animation!.wiggle({ mode: 'dynamics', amplitude: 2 }));
    expect(v.managers.history.steps.value.length).toBe(steps + 1);

    expectOk(await v.api.history.undo());
    expect(wiggleOf(v.plugin)).toMatchObject({ mode: 'dynamics', amplitude: 0.5 });
  });

  it('wiggles by uncertainty with layers, and clears them', async () => {
    const v = await headlessViewer({ files: [cbs()], extensions: [AnimationFeature] });
    expectOk(await v.api.animation!.wiggle({ mode: 'uncertainty' }));
    expect(wiggleOf(v.plugin).mode).toBe('uncertainty');
    expect(wiggleCells(v.plugin).length).toBeGreaterThan(0);

    expectOk(await v.api.animation!.wiggle({ mode: 'off' }));
    expect(wiggleOf(v.plugin).mode).toBe('off');
    expect(wiggleCells(v.plugin).length).toBe(0);
    expect(await v.api.animation!.wiggle({ mode: 'off' })).toMatchObject({ ok: true, unchanged: true });

    expectOk(await v.api.history.undo());
    expect(wiggleOf(v.plugin).mode).toBe('uncertainty');
  });

  it('says when the uncertainty wiggle is already so, and what it did in the past tense', async () => {
    const v = await headlessViewer({ files: [cbs()], extensions: [AnimationFeature] });
    const first = await v.api.animation!.wiggle({ mode: 'uncertainty', amplitude: 1 });
    expectOk(first);
    expect(first.did).toBe('Started the wiggle (uncertainty, amplitude 1).');
    const steps = v.managers.history.steps.value.length;
    // A repeat does not make another step.
    expect(await v.api.animation!.wiggle({ mode: 'uncertainty', amplitude: 1 })).toMatchObject({
      ok: true,
      unchanged: true,
    });
    expect(v.managers.history.steps.value.length).toBe(steps);
    // Another amplitude is another wiggle.
    expectOk(await v.api.animation!.wiggle({ mode: 'uncertainty', amplitude: 2 }));
    expect(v.managers.history.steps.value.length).toBe(steps + 1);
    expect(await v.api.animation!.wiggle({ mode: 'uncertainty', amplitude: 2 })).toMatchObject({ unchanged: true });

    const dynamics = await v.api.animation!.wiggle({ mode: 'dynamics', amplitude: 1 });
    expectOk(dynamics);
    expect(dynamics.did).toBe('Started the wiggle (dynamics, amplitude 1).');
    const off = await v.api.animation!.wiggle({ mode: 'off' });
    expectOk(off);
    expect(off.did).toBe('Stopped the wiggle.');
    // One name in the command menu and the History.
    expect(v.api.animation!.wiggle.definition.title).toBe('Wiggle');
    expect(v.managers.commands.get('animation.wiggle')?.title).toBe('Wiggle…');
  });

  it('switches from one mode to the other', async () => {
    const v = await headlessViewer({ files: [cbs()], extensions: [AnimationFeature] });
    expectOk(await v.api.animation!.wiggle({ mode: 'uncertainty' }));
    expectOk(await v.api.animation!.wiggle({ mode: 'dynamics' }));
    expect(wiggleOf(v.plugin).mode).toBe('dynamics');
    expect(wiggleCells(v.plugin).length).toBe(0);
  });

  it('sets the type, speed, frequency and tumble; keeps what is not given', async () => {
    const v = await headlessViewer({ files: [cbs()], extensions: [AnimationFeature] });
    expectOk(
      await v.api.animation!.wiggle({ mode: 'dynamics', amplitude: 0.5, type: 'group', speed: 3, frequency: 0.5 }),
    );
    expect(wiggleOf(v.plugin)).toMatchObject({
      mode: 'dynamics',
      amplitude: 0.5,
      type: 'group',
      speed: 3,
      frequency: 0.5,
    });
    const repr = () =>
      (
        v.plugin.managers.structure.hierarchy.current.structures[0].components[0].representations[0].cell.transform
          .params as any
      ).type.params.animation;
    expect(repr()).toMatchObject({ wiggleMode: 'group', wiggleSpeed: 3, wiggleFrequency: 0.5 });

    // the tumble alone: the wiggle as it was
    const t = await v.api.animation!.wiggle({ tumble: 2, tumbleSpeed: 4 });
    expectOk(t);
    expect(t.did).toBe('Changed the wiggle (dynamics, amplitude 0.5); the copies tumble (2 Å).');
    expect(wiggleOf(v.plugin)).toMatchObject({ mode: 'dynamics', amplitude: 0.5, speed: 3, tumble: 2, tumbleSpeed: 4 });
    expect(repr()).toMatchObject({ tumbleAmplitude: 2, tumbleSpeed: 4 });

    // off stops the wiggle, the tumble goes on; then it stops too
    expect(await v.api.animation!.wiggle({ mode: 'off' })).toMatchObject({
      did: 'Stopped the wiggle; the copies tumble (2 Å).',
    });
    expect(await v.api.animation!.wiggle({ mode: 'off' })).toMatchObject({ unchanged: true });
    expect(await v.api.animation!.wiggle({ tumble: 0 })).toMatchObject({ did: 'Stopped the tumble.' });
    expect(wiggleOf(v.plugin)).toMatchObject({ mode: 'off', tumble: 0 });
  });

  it('needs a structure', async () => {
    const v = await headlessViewer({ extensions: [AnimationFeature] });
    expect(await v.api.animation!.wiggle({ mode: 'dynamics' })).toMatchObject({ ok: false, code: 'not-found' });
  });
});
