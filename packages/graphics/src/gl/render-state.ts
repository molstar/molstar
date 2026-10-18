/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Backend-neutral fixed-function state, named after WebGPU. WebGL applies it through the
 * cached `WebGLState` setters, WebGPU makes it part of the pipeline key.
 */

export type BlendFactor =
  | 'zero'
  | 'one'
  | 'src'
  | 'one-minus-src'
  | 'src-alpha'
  | 'one-minus-src-alpha'
  | 'dst'
  | 'one-minus-dst'
  | 'dst-alpha'
  | 'one-minus-dst-alpha'
  /** the blend constant, see `GraphicsContext.setBlendConstant` */
  | 'constant'
  | 'one-minus-constant';

/** `min` and `max` ignore the factors, use `one` for both (required by WebGPU) */
export type BlendOperation = 'add' | 'subtract' | 'reverse-subtract' | 'min' | 'max';

export interface BlendComponent {
  readonly operation: BlendOperation;
  readonly srcFactor: BlendFactor;
  readonly dstFactor: BlendFactor;
}

export interface BlendState {
  readonly color: BlendComponent;
  readonly alpha: BlendComponent;
}

export type CompareFunction =
  | 'never'
  | 'less'
  | 'equal'
  | 'less-equal'
  | 'greater'
  | 'not-equal'
  | 'greater-equal'
  | 'always';

export interface DepthState {
  readonly compare: CompareFunction;
  readonly write: boolean;
}

export type StencilOperation =
  | 'keep'
  | 'zero'
  | 'replace'
  | 'invert'
  | 'increment-clamp'
  | 'decrement-clamp'
  | 'increment-wrap'
  | 'decrement-wrap';

export interface StencilFaceState {
  readonly compare: CompareFunction;
  readonly failOp: StencilOperation;
  readonly depthFailOp: StencilOperation;
  readonly passOp: StencilOperation;
}

export interface StencilState {
  readonly front: StencilFaceState;
  readonly back: StencilFaceState;
  readonly readMask: number;
  readonly writeMask: number;
  /** dynamic state on WebGPU */
  readonly reference: number;
}

export type CullMode = 'none' | 'front' | 'back';
export type FrontFace = 'ccw' | 'cw';
export type ColorWriteMask = readonly [red: boolean, green: boolean, blue: boolean, alpha: boolean];

/** Unset fields use the defaults: no blending, all channels written, no depth or stencil test, no culling, `ccw` */
export interface RenderState {
  readonly blend?: BlendState;
  readonly colorWrite?: ColorWriteMask;
  /** depth test and write, both off when unset */
  readonly depth?: DepthState;
  readonly stencil?: StencilState;
  readonly cull?: CullMode;
  readonly frontFace?: FrontFace;
}

/** Values to clear the attachments of the bound target to, within the scissor rectangle and ignoring write masks */
export interface ClearValues {
  /** integer attachments are cleared to the truncated values */
  readonly color?: readonly [red: number, green: number, blue: number, alpha: number];
  readonly depth?: number;
  readonly stencil?: number;
}

function BlendComponent(operation: BlendOperation, srcFactor: BlendFactor, dstFactor: BlendFactor): BlendComponent {
  return { operation, srcFactor, dstFactor };
}

export function BlendState(color: BlendComponent, alpha: BlendComponent = color): BlendState {
  return { color, alpha };
}

export namespace BlendState {
  export const Component = BlendComponent;

  /** `src * 1 + dst * 1` */
  export const Additive = BlendState(BlendComponent('add', 'one', 'one'));
  /** `max(src, dst)` */
  export const Max = BlendState(BlendComponent('max', 'one', 'one'));
  /** `src * 1 + dst * (1 - srcAlpha)`, for premultiplied colors */
  export const Premultiplied = BlendState(BlendComponent('add', 'one', 'one-minus-src-alpha'));
  /** straight alpha for color, `src * 1 + dst * (1 - srcAlpha)` for alpha */
  export const Alpha = BlendState(
    BlendComponent('add', 'src-alpha', 'one-minus-src-alpha'),
    BlendComponent('add', 'one', 'one-minus-src-alpha'),
  );
  /** `src * constant + dst * (1 - constant)` */
  export const Constant = BlendState(BlendComponent('add', 'constant', 'one-minus-constant'));
}

export namespace RenderState {
  /** Fullscreen and compute passes: no blending, depth, stencil or culling, all channels written */
  export const Default: RenderState = {};
  export const Additive: RenderState = { blend: BlendState.Additive };
  export const Premultiplied: RenderState = { blend: BlendState.Premultiplied };
  export const Alpha: RenderState = { blend: BlendState.Alpha };
  export const Constant: RenderState = { blend: BlendState.Constant };
}

export interface TargetRect {
  x: number;
  y: number;
  width: number;
  height: number;
}

/**
 * The part of a scissor rectangle (bottom-left origin) inside a target of `targetWidth` x `targetHeight`, empty where
 * they do not overlap: WebGPU requires scissors within the target, GL clips them to it. Viewports need no clamping,
 * both backends crop a viewport partly outside the target
 */
export function intersectTargetRect(
  out: TargetRect,
  x: number,
  y: number,
  width: number,
  height: number,
  targetWidth: number,
  targetHeight: number,
): TargetRect {
  const x0 = Math.min(Math.max(x, 0), targetWidth);
  const y0 = Math.min(Math.max(y, 0), targetHeight);
  out.x = x0;
  out.y = y0;
  out.width = Math.max(0, Math.min(x + width, targetWidth) - x0);
  out.height = Math.max(0, Math.min(y + height, targetHeight) - y0);
  return out;
}
