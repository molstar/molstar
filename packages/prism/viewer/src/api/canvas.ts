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
import type { ReadonlyState } from '@molstar/prism-ui/core';
import type { CanvasState } from '../sources/viewer-sources.js';
import { type DataPalette, dataPaletteOf, setDataPalette } from '../sources/palette.js';
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
      if (ctx.sources.canvas.value?.background === hex) return ApiResult.already(`The background is already ${color}.`);
      ctx.plugin.setCanvas3DProps({ renderer: { backgroundColor: value } });
      return ApiResult.ok(`Background set to ${color}.`, { view: true });
    },
  });

  // README §3.7: the categorical themes' palette, colorblind-safe or not; the representations on the old default follow.
  readonly palette = this.operation({
    id: 'canvas.palette',
    title: 'Set the data palette',
    params: {
      palette: PD.Select<DataPalette>(
        'colorblind',
        [
          ['default', 'Default'],
          ['colorblind', 'Colorblind-safe (Okabe–Ito)'],
        ],
        { description: 'The colors of chains, entities, polymers and other categories' },
      ),
    },
    async run({ palette }, ctx) {
      if (dataPaletteOf(ctx.plugin) === palette) return ApiResult.already(`The data palette is already ${palette}.`);
      const nodes = (await setDataPalette(ctx.plugin, palette)) ?? [];
      return ApiResult.ok(
        palette === 'colorblind'
          ? 'Chains and entities now use colorblind-safe colors.'
          : 'Chains and entities use the default colors.',
        { nodes },
      );
    },
  });

  readonly get = this.operation({
    id: 'canvas.get',
    title: 'Describe the canvas',
    params: {},
    read: true,
    async run(_, ctx) {
      const canvas = ctx.sources.canvas.value;
      const p = ctx.plugin.canvas3d?.props;
      if (!canvas || !p) return ApiResult.ok('No canvas (a headless viewer).', {}, { value: undefined });
      // What the canvas operations set, so a script reads it back as it would give it (operations.md §1).
      const on = (m: { name: string }) => m.name === 'on';
      const value = {
        background: canvas.background,
        transparent: canvas.transparentBackground,
        illumination: canvas.illumination,
        lighting: {
          occlusion: on(p.postprocessing.occlusion),
          shadow: on(p.postprocessing.shadow),
          outline: on(p.postprocessing.outline),
          fog: on(p.cameraFog),
          illumination: p.illumination.enabled,
        },
        camera: { mode: p.camera.mode, fov: p.camera.fov },
        clipping: { clip: p.cameraClipping.radius, far: p.cameraClipping.far },
        xr: !!ctx.plugin.canvas3d?.xr?.isPresenting.value,
      };
      const lit = (Object.keys(value.lighting) as (keyof typeof value.lighting)[]).filter((k) => value.lighting[k]);
      return ApiResult.ok(
        `Background ${canvas.background}${canvas.transparentBackground ? ' (transparent)' : ''}; lighting: ${lit.length ? lit.join(', ') : 'plain'}; ${p.camera.mode} camera, ${p.camera.fov}°${p.cameraClipping.radius ? `; ${p.cameraClipping.radius}% clipped` : ''}.`,
        {},
        { value },
      );
    },
  });
}
