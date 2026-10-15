/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import { decodeColor } from '@molstar/core/util/color/utils';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ReadonlyState } from '@molstar/prism-ui';
import type { CanvasState } from '../sources/viewer-sources.js';
import { ApiArea, ApiResult } from './area.js';

// `canvas` (operations.md §2): background now; illumination, lighting, clipping and XR with `canvas-settings` (plan
// phase 3). Canvas props are set through the plugin's emitting path, so the canvas source follows (state.md §3).

export class CanvasApi extends ApiArea {
  readonly id = 'canvas';

  /** The canvas as the viewer reads it; undefined without one (headless). */
  get state(): ReadonlyState<CanvasState | undefined> {
    return this.ctx.sources.canvas;
  }

  readonly background = this.operation({
    id: 'canvas.background',
    title: 'Set the background',
    params: { color: PD.Text('#ffffff', { description: 'A CSS color: white, #ffffff, rgb(255, 255, 255)' }) },
    modules: ['canvas'],
    async run({ color }, ctx) {
      const value = decodeColor(color.trim());
      if (value === undefined) {
        return ApiResult.fail(
          'invalid',
          `'${color}' is not a color.`,
          'Use a CSS color: white, #ffffff, rgb(255, 255, 255).',
        );
      }
      const hex = Color.toHexStyle(value);
      if (ctx.sources.canvas.value?.background === hex) return ApiResult.ok(`The background is already ${color}.`);
      ctx.plugin.setCanvas3DProps({ renderer: { backgroundColor: value } });
      return ApiResult.ok(`Background set to ${color}.`, { view: true });
    },
  });

  readonly get = this.operation({
    id: 'canvas.get',
    title: 'Describe the canvas',
    params: {},
    read: true,
    async run(_, ctx) {
      const canvas = ctx.sources.canvas.value;
      if (!canvas) return ApiResult.ok('No canvas (a headless viewer).', {}, { value: undefined });
      const value = {
        background: canvas.background,
        transparent: canvas.transparentBackground,
        illumination: canvas.illumination,
      };
      return ApiResult.ok(
        `Background ${canvas.background}${canvas.transparentBackground ? ' (transparent)' : ''}; illumination ${canvas.illumination ? 'on' : 'off'}.`,
        {},
        { value },
      );
    },
  });
}
