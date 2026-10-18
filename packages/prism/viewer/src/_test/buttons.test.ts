/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { headlessViewer } from './harness.js';

// `options.ui.buttons` (viewer-spec.md §2.2): each option hides the controls it names, by their ids in the shell.

describe('options.ui.buttons', () => {
  it('hides the toolbar items and canvas controls they name', async () => {
    const off = {
      reset: false,
      screenshot: false,
      expand: false,
      fullscreen: false,
      settings: false,
      selectionMode: false,
      animation: false,
      trajectory: false,
    };
    const v = await headlessViewer({ spec: { options: { ui: { buttons: off } } } });
    expect(v.managers.viewState.layout.value.buttons).toEqual({
      reset: false,
      export: false,
      expand: false,
      'canvas-settings.toolbar.canvas': false,
      'selection-mode': false,
      'animation.toolbar.animation': false,
      trajectory: false,
    });
  });
});
