/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Frame } from '@molstar/graphics/gl/renderable';
import type { Renderer } from '@molstar/graphics/gl/renderer';
import type { Scene } from '@molstar/graphics/gl/scene';
import type { Camera } from './camera.js';
import type { StereoCamera } from './camera/stereo.js';
import type { Canvas3DAttribs, Canvas3DProps } from './canvas3d.js';
import type { Helper } from './helper/helper.js';
import type { PickHelper } from './helper/pick-helper.js';
import type { ShaderManager } from './helper/shader-manager.js';
import type { HiZPass } from './passes/hi-z.js';
import type { MultiSampleHelper } from './passes/multi-sample.js';
import type { Passes } from './passes/passes.js';

/**
 * The parts of a canvas (`Canvas3D.create`) its modules work with, read when used: a context change (`setContext`)
 * replaces the graphics context and what is bound to it (the renderer, passes and helpers), the frame advances
 */
export interface Canvas3DParts {
  readonly graphics: GraphicsContext;
  readonly renderer: Renderer;
  readonly passes: Passes;
  readonly hiZ: HiZPass;
  readonly multiSampleHelper: MultiSampleHelper;
  readonly pickHelper: PickHelper;
  readonly shaderManager: ShaderManager;
  readonly scene: Scene;
  readonly helper: Helper;
  readonly camera: Camera;
  readonly stereoCamera: StereoCamera;
  /** the frame being drawn (`Renderer.update`) */
  readonly frame: Frame;
  /** the canvas's props, changed in place by `setProps` */
  readonly props: Canvas3DProps;
  readonly attribs: Canvas3DAttribs;
  /** `dispose` was called */
  readonly disposed: boolean;
  /** the context is lost (nothing to plan or compile until it is restored or replaced) */
  readonly contextLost: boolean;
}
