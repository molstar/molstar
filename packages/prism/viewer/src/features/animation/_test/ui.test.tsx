/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { afterEach, describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from '../../../_test/harness.js';
import { AnimationFeature } from '../index.js';
import { animationUI } from '../ui.js';

// The Animation card's sections over a headless viewer: the part ids, what they render, and that the wiggle acts through
// the API.

afterEach(() => {
  document.body.innerHTML = '';
});

describe('the Animation card', () => {
  it('has the sections, the overlay and the toolbar item under the feature id', async () => {
    const v = await headlessViewer({ files: [cbs()], extensions: [AnimationFeature] });
    const ui = animationUI({ api: v.api.as({ kind: 'person' }), sources: v.sources, managers: v.managers });
    expect(ui.sections?.map((s) => s.id)).toEqual(['animation.section.play', 'animation.section.wiggle']);
    expect(ui.sections?.every((s) => s.slot === 'animation')).toBe(true);
    expect(ui.overlays?.map((o) => o.id)).toEqual(['animation.overlay.animation']);
    expect(ui.toolbar?.[0]).toMatchObject({
      id: 'animation.toolbar.animation',
      icon: 'clapperboard',
      group: 'scene',
      overlay: 'animation.overlay.animation',
    });

    const mount = (i: number) => {
      const host = document.createElement('div');
      document.body.append(host);
      ui.sections![i].render.mount(host, { shell: new ShellModel() } as never, undefined as never);
      return host;
    };
    const play = mount(0);
    expect(play.textContent).toContain('Animation');
    expect(play.textContent).toContain('Camera Spin');
    const wiggle = mount(1);
    const apply = [...wiggle.querySelectorAll('button')].find((b) => b.textContent?.includes('Apply'))!;
    apply.click();
    await until(() => v.api.animation!.state.value.wiggle.mode === 'dynamics');
    await until(() => wiggle.textContent?.includes('Wiggling by dynamics') ?? false);
    const clear = [...wiggle.querySelectorAll('button')].find((b) => b.textContent?.includes('Clear'))!;
    clear.click();
    await until(() => v.api.animation!.state.value.wiggle.mode === 'off');
    await until(() => v.managers.history.running.value === 0);
    expectOk(await v.api.history.undo());
  });
});
