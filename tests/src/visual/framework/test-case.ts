/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Canvas3D, Canvas3DContext, type PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { type GraphicsContext } from '@molstar/graphics/gl/context';
import { type WebGLContext } from '@molstar/graphics/gl/webgl/context';
import { Representation } from '@molstar/graphics/repr/representation';
import { type CapabilityKey } from './capabilities.js';
import type { BackendProfile } from './backend.js';
import { type Tolerance } from './compare.js';

export interface SceneContext {
  /** The graphics context of the backend under test */
  readonly graphics: GraphicsContext;
  /** The WebGL context, for cases that build GPU data with WebGL; such cases are skipped on other backends */
  readonly webgl: WebGLContext;
  readonly canvas3d: Canvas3D;
}

export interface TestView {
  /** Viewing direction, from camera towards the scene center; default [0, 0, -1] */
  direction?: [number, number, number];
  up?: [number, number, number];
  /** Values > 1 zoom in; default 1 */
  zoom?: number;
}

export interface VisualTestCase {
  /** Stable, unique, `group/sub/name`; also used for reference file names */
  readonly id: string;
  readonly description?: string;
  readonly tags?: readonly string[];
  readonly requires?: readonly CapabilityKey[];
  /**
   * Documented fallback for profiles without `fallback.requires`: such a profile renders the case as is
   * (the feature falls back internally), and is compared against a profile with the capability rendering
   * the case with `fallback.overrides` applied
   */
  readonly fallback?: {
    readonly requires: readonly CapabilityKey[];
    readonly overrides: Partial<Pick<VisualTestCase, 'contextProps' | 'canvasProps' | 'scene'>>;
  };
  readonly contextProps?: Partial<Pick<Canvas3DContext.Props, 'transparency'>>;
  readonly canvasProps?: PartialCanvas3DProps;
  /** Camera is fitted to the scene bounding sphere, then `camera` overrides are applied */
  readonly view?: TestView;
  readonly camera?: Partial<Camera.Snapshot>;
  readonly size?: { width: number; height: number };
  readonly tolerance?: Partial<Tolerance>;
  /** Render with a fresh backend instance instead of the shared one */
  readonly isolate?: boolean;
  /** Canvas3D frames drawn before the capture, for passes that only run interactively (e.g. hi-z) */
  readonly frames?: number;
  /**
   * The context released before the capture (`Canvas3DContext.release`), checked to stay lost, then restored
   * (`restore`): the capture is of the scene built again on the new context; needs `releaseContext`
   */
  readonly releaseContext?: boolean;
  /**
   * `image-pass` (default): an image pass renders the scene at the full size; `canvas`: a Canvas3D frame read from its
   * drawing buffer, for what only Canvas3D draws (viewports other than the canvas, stereo, Hi-Z culling)
   */
  readonly capture?: 'image-pass' | 'canvas';
  /** Called once per backend; GPU resources are per context */
  readonly scene: (ctx: SceneContext) => Promise<Representation.Any[]>;
  /**
   * Called with the scene added and the camera fitted, before the capture (e.g. picking with `canvas3d`, marking what
   * was picked): the errors it returns fail the case
   */
  readonly interact?: (ctx: SceneContext, reprs: Representation.Any[]) => Promise<string[]>;
  /**
   * A check without an image (e.g. compute kernels against CPU results): run on the target backend only, the case
   * passes if it returns no errors; `requires` applies to the target only and `scene` is not rendered; `profile` makes
   * contexts of its own (`createContext`)
   */
  readonly check?: (graphics: GraphicsContext, profile: BackendProfile) => Promise<string[]>;
}

export const DefaultTestSize = { width: 256, height: 256 };

export function testIdToFileName(id: string) {
  return `${id.replace(/\//g, '.')}.png`;
}

/** `*` matches within a path segment, `**` across segments */
export function globToRegExp(glob: string) {
  const escaped = glob.replace(/[.+^${}()|[\]\\]/g, '\\$&');
  const pattern = escaped
    .replace(/\*\*/g, '\u0000')
    .replace(/\*/g, '[^/]*')
    .replace(/\u0000/g, '.*');
  return new RegExp(`^${pattern}$`);
}

export function matchesAny(id: string, globs: readonly string[]) {
  return globs.some((g) => globToRegExp(g).test(id));
}
