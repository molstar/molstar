/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createRenderObject } from '../render-object.js';
import { Scene } from '../scene.js';
import { GraphicsRenderVariants } from '../render-item.js';
import { getGLContext, tryGetGLContext } from './gl.js';
import { setDebugMode } from '@molstar/core/util/debug';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Clip } from '@molstar/core/util/clip';
import { DirectVolume } from '@molstar/graphics/geo/geometry/direct-volume/direct-volume';

export function createDirectVolume(
  renderMode: PD.Values<DirectVolume.Params>['renderMode'] = 'composite',
  clipped = false,
) {
  const directVolume = DirectVolume.createEmpty();
  const props = { ...PD.getDefaultValues(DirectVolume.Params), renderMode };
  // a clip plane: the pixel clipping of the ray marcher
  if (clipped) props.clip = { variant: 'pixel', objects: [Clip.Params.objects.ctor()] };
  const values = DirectVolume.Utils.createValuesSimple(directVolume, props, ColorNames.orange, 1);
  const state = DirectVolume.Utils.createRenderableState(props);
  return createRenderObject('direct-volume', values, state, -1);
}

describe('direct-volume', () => {
  const ctx = tryGetGLContext(32, 32);

  for (const renderMode of ['composite', 'maximum', 'attenuated', 'additive'] as const) {
    for (const clipped of [false, true]) {
      (ctx ? it : it.skip)(`basic (${renderMode}${clipped ? ', clipped' : ''})`, { tags: ['medium'] }, async () => {
        const ctx = getGLContext(32, 32);
        const scene = Scene.create(ctx);
        const directVolume = createDirectVolume(renderMode, clipped);
        scene.add(directVolume);
        setDebugMode(true);
        expect(() => scene.commit()).not.toThrow();
        // programs link on first use: link each variant, so that a shader error throws (debug mode)
        for (const r of scene.renderables) for (const v of GraphicsRenderVariants) r.getProgram(v).finalize(true);
        setDebugMode(false);
        ctx.destroy();
      });
    }
  }
});
