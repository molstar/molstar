/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ensemble } from '../../../api/_test/trajectory-fixtures.js';
import { cbs, expectOk, headlessViewer, until } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import { AnimationFeature } from '../index.js';
import { createAnimationModel } from '../model.js';

// The Animation card's model: the choice of animation with its params, why it cannot play, and the wiggle's form.

describe('the animation model', () => {
  it('forms the animations as one choice with their params, and says why one cannot play', async () => {
    const v = await headlessViewer({ files: [cbs()], extensions: [AnimationFeature] });
    const model = createAnimationModel(v.api.as({ kind: 'person' }), v.api.animation!);
    const form = model.form!;
    const choice = form.params.animation as any;
    expect(choice.select.options.map((o: string[]) => o[0])).toContain('built-in.animate-camera-spin');
    // Each option's params are the animation's own.
    expect(Object.keys(choice.map('built-in.animate-camera-spin').params)).toEqual(['durationInMs', 'speed', 'axis']);

    form.set('animation', {
      name: 'built-in.animate-model-index',
      params: choice.map('built-in.animate-model-index').defaultValue,
    } as never);
    expect(model.blocked.value).toMatch(/No trajectory/);
    form.set('animation', {
      name: 'built-in.animate-camera-spin',
      params: choice.map('built-in.animate-camera-spin').defaultValue,
    } as never);
    expect(model.blocked.value).toBeUndefined();
  });

  it('builds the form again when a structure to unwind loads, keeping the choice and the values that still apply', async () => {
    const v = await headlessViewer({ extensions: [AnimationFeature] });
    const model = createAnimationModel(v.api.as({ kind: 'person' }), v.api.animation!);
    const unwind = 'built-in.animate-assembly-unwind';
    const targets = () =>
      ((model.form!.params.animation as any).map(unwind).params.target.options as string[][]).map((o) => o[1]);
    expect(targets()).toEqual(['All']);
    const defaults = (model.form!.params.animation as any).map(unwind).defaultValue;
    model.form!.set('animation', { name: unwind, params: { ...defaults, durationInMs: 1200 } } as never);
    const before = model.form;
    // Nothing changed: the same form.
    model.refresh();
    expect(model.form).toBe(before);

    expectOk(await v.api.load.files({ files: [cbs()] }));
    model.refresh();
    expect(model.form).not.toBe(before);
    expect(targets()).toHaveLength(2);
    expect(model.form!.values.value.animation).toMatchObject({
      name: unwind,
      params: { durationInMs: 1200, target: 'all' },
    });
  });

  it('plays and stops through the API, and applies and clears the wiggle', async () => {
    const v = await createPrismViewerHeadless({
      spec: {
        ...DefaultPrismViewerSpec(),
        // The default spec has the feature (in the place it lists it).
        extensions: DefaultPrismViewerSpec().extensions ?? [],
      },
      modules: ['canvas', 'files'],
    });
    try {
      expectOk(await v.api.load.files({ files: [ensemble()] }));
      const model = createAnimationModel(v.api.as({ kind: 'person' }), v.api.animation!);
      const choice = model.form!.params.animation as any;
      model.form!.set('animation', {
        name: 'built-in.animate-model-index',
        params: choice.map('built-in.animate-model-index').defaultValue,
      } as never);
      expect(model.blocked.value).toBeUndefined();

      model.play();
      await until(() => model.isPlaying.value);
      model.stop();
      await until(() => !model.isPlaying.value);

      model.wiggleForm.set('amplitude', 2);
      model.wiggleForm.set('speed', 4);
      model.wiggleForm.set('tumble', { amplitude: 1.5, speed: 2, frequency: 0.2 });
      model.applyWiggle();
      await until(() => model.wiggle.value.mode === 'dynamics');
      expect(model.wiggle.value).toMatchObject({ amplitude: 2, speed: 4, tumble: 1.5, tumbleSpeed: 2 });
      model.clearWiggle();
      await until(() => model.wiggle.value.mode === 'off' && model.wiggle.value.tumble === 0);
      // set elsewhere: the form takes it when the section shows again
      expectOk(await v.api.animation!.wiggle({ mode: 'dynamics', amplitude: 0.7, type: 'group' }));
      model.refreshWiggle();
      expect(model.wiggleForm.values.value).toMatchObject({ mode: 'dynamics', amplitude: 0.7, type: 'group' });
    } finally {
      await v.dispose();
    }
  });

  it('labels the three numbers of a vector, and plays with them as a vector', async () => {
    const v = await createPrismViewerHeadless({
      spec: { ...DefaultPrismViewerSpec(), extensions: DefaultPrismViewerSpec().extensions ?? [] },
      modules: ['canvas', 'files'],
    });
    try {
      const model = createAnimationModel(v.api.as({ kind: 'person' }), v.api.animation!);
      const choice = model.form!.params.animation as any;
      const params = choice.map('built-in.animate-camera-spin').params;
      expect(params.durationInMs.label).toBe('Duration (ms)');
      expect(params.axis.type).toBe('group');
      expect(Object.values(params.axis.params).map((p: any) => p.label)).toEqual(['Axis x', 'Axis y', 'Axis z']);

      const defaults = choice.map('built-in.animate-camera-spin').defaultValue;
      model.form!.set('animation', {
        name: 'built-in.animate-camera-spin',
        params: { ...defaults, axis: { x: 0, y: 1, z: 0 } },
      } as never);
      model.play();
      await until(() => model.isPlaying.value);
      expect(v.plugin.managers.animation.current?.paramValues.axis).toEqual([0, 1, 0]);
      model.stop();
      await until(() => !model.isPlaying.value);
    } finally {
      await v.dispose();
    }
  });
});
