/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Renderer, type RendererProps } from '@molstar/graphics/gl/renderer';
import type { OcclusionPyramid } from '@molstar/graphics/gl/compute/cell-cull';
import type { Scene } from '@molstar/graphics/gl/scene';
import type { Viewport } from './camera/util.js';
import type { Helper } from './helper/helper.js';
import { PickHelper } from './helper/pick-helper.js';
import { RayHelper } from './helper/ray-helper.js';
import { ShaderManager } from './helper/shader-manager.js';
import { HiZPass, type HiZProps } from './passes/hi-z.js';
import { MultiSampleHelper } from './passes/multi-sample.js';
import type { Passes } from './passes/passes.js';
import type { PickOptions } from './passes/pick.js';

/** The parts of a canvas bound to its graphics context, made with it and again by a context change */
export interface Canvas3DContextParts {
  readonly hiZ: HiZPass;
  readonly renderer: Renderer;
  readonly shaderManager: ShaderManager;
  readonly pickHelper: PickHelper;
  readonly rayHelper: RayHelper;
  readonly multiSampleHelper: MultiSampleHelper;
}

/** What the context's parts are made of: the context and passes, and what of the canvas they keep */
export interface Canvas3DContextPartsInputs {
  readonly graphics: GraphicsContext;
  readonly passes: Passes;
  readonly canvas: HTMLCanvasElement | undefined;
  readonly scene: Scene;
  readonly helper: Helper;
  readonly hiZProps: Partial<HiZProps>;
  readonly rendererProps: Partial<RendererProps>;
  readonly pickOptions: PickOptions;
  readonly viewport: Viewport;
  /** the depth pyramid of the previous frame, for GPU culling (`Renderer.setOcclusionPyramid`) */
  readonly getOcclusionPyramid: () => OcclusionPyramid | undefined;
}

/** Makes the parts of a canvas bound to `inputs.graphics`, the renderer culling by Hi-Z's occlusion */
export function createCanvas3DContextParts(inputs: Canvas3DContextPartsInputs): Canvas3DContextParts {
  const { graphics, passes, canvas, scene, helper, pickOptions, viewport } = inputs;
  const hiZ = new HiZPass(graphics, passes.draw, canvas, inputs.hiZProps);
  const renderer = Renderer.create(graphics, inputs.rendererProps);
  renderer.setOcclusionTest(hiZ.occlusionTest);
  renderer.setOcclusionPyramid(inputs.getOcclusionPyramid);
  return {
    hiZ,
    renderer,
    shaderManager: new ShaderManager(graphics, scene, helper.overlay),
    pickHelper: new PickHelper(graphics, renderer, scene, helper, passes.pick, viewport, pickOptions),
    rayHelper: new RayHelper(graphics, renderer, scene, helper, pickOptions),
    multiSampleHelper: new MultiSampleHelper(passes.multiSample),
  };
}

/** Releases what the context's parts hold of their context (before it is replaced) */
export function disposeCanvas3DContextParts(parts: Canvas3DContextParts) {
  parts.renderer.dispose();
  parts.hiZ.dispose();
  parts.pickHelper.dispose();
  parts.rayHelper.dispose();
}
