/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type {
  BlendFactor,
  BlendOperation,
  ClearValues,
  CompareFunction,
  RenderState,
  StencilOperation,
} from '../render-state.js';
import { type GLRenderingContext, isWebGL2 } from './compat.js';
import type { WebGLExtensions } from './extensions.js';
import { hasIntegerColorAttachment } from './framebuffer.js';
import type { WebGLState } from './state.js';
import { checkRenderPass } from './render-pass.js';

function getBlendFactor(gl: GLRenderingContext, factor: BlendFactor) {
  switch (factor) {
    case 'zero':
      return gl.ZERO;
    case 'one':
      return gl.ONE;
    case 'src':
      return gl.SRC_COLOR;
    case 'one-minus-src':
      return gl.ONE_MINUS_SRC_COLOR;
    case 'src-alpha':
      return gl.SRC_ALPHA;
    case 'one-minus-src-alpha':
      return gl.ONE_MINUS_SRC_ALPHA;
    case 'dst':
      return gl.DST_COLOR;
    case 'one-minus-dst':
      return gl.ONE_MINUS_DST_COLOR;
    case 'dst-alpha':
      return gl.DST_ALPHA;
    case 'one-minus-dst-alpha':
      return gl.ONE_MINUS_DST_ALPHA;
    case 'constant':
      return gl.CONSTANT_COLOR;
    case 'one-minus-constant':
      return gl.ONE_MINUS_CONSTANT_COLOR;
  }
}

function getBlendOperation(gl: GLRenderingContext, extensions: WebGLExtensions, operation: BlendOperation) {
  switch (operation) {
    case 'add':
      return gl.FUNC_ADD;
    case 'subtract':
      return gl.FUNC_SUBTRACT;
    case 'reverse-subtract':
      return gl.FUNC_REVERSE_SUBTRACT;
    case 'min':
      return extensions.blendMinMax.MIN;
    case 'max':
      return extensions.blendMinMax.MAX;
  }
}

function getCompareFunction(gl: GLRenderingContext, compare: CompareFunction) {
  switch (compare) {
    case 'never':
      return gl.NEVER;
    case 'less':
      return gl.LESS;
    case 'equal':
      return gl.EQUAL;
    case 'less-equal':
      return gl.LEQUAL;
    case 'greater':
      return gl.GREATER;
    case 'not-equal':
      return gl.NOTEQUAL;
    case 'greater-equal':
      return gl.GEQUAL;
    case 'always':
      return gl.ALWAYS;
  }
}

function getStencilOperation(gl: GLRenderingContext, operation: StencilOperation) {
  switch (operation) {
    case 'keep':
      return gl.KEEP;
    case 'zero':
      return gl.ZERO;
    case 'replace':
      return gl.REPLACE;
    case 'invert':
      return gl.INVERT;
    case 'increment-clamp':
      return gl.INCR;
    case 'decrement-clamp':
      return gl.DECR;
    case 'increment-wrap':
      return gl.INCR_WRAP;
    case 'decrement-wrap':
      return gl.DECR_WRAP;
  }
}

export interface WebGLRenderStateApplier {
  setRenderState: (renderState: RenderState) => void;
  /** Masks depth and stencil writes (a pass with `depthReadOnly`), until set to false */
  setDepthReadOnly: (readOnly: boolean) => void;
  setViewport: (x: number, y: number, width: number, height: number) => void;
  setScissor: (x: number, y: number, width: number, height: number) => void;
  disableScissor: () => void;
  clear: (values: ClearValues) => void;
}

export function createRenderStateApplier(
  gl: GLRenderingContext,
  state: WebGLState,
  extensions: WebGLExtensions,
): WebGLRenderStateApplier {
  const intClearValue = new Int32Array(4);
  /** to restore the write masks after clearing */
  let current: RenderState = {};
  let depthReadOnly = false;

  return {
    setDepthReadOnly: (readOnly: boolean) => {
      depthReadOnly = readOnly;
    },
    setRenderState: (renderState: RenderState) => {
      current = renderState;
      const { blend, colorWrite, depth, stencil, cull, frontFace } = renderState;

      if (blend) {
        state.enable(gl.BLEND);
        const { color, alpha } = blend;
        state.blendFuncSeparate(
          getBlendFactor(gl, color.srcFactor),
          getBlendFactor(gl, color.dstFactor),
          getBlendFactor(gl, alpha.srcFactor),
          getBlendFactor(gl, alpha.dstFactor),
        );
        state.blendEquationSeparate(
          getBlendOperation(gl, extensions, color.operation),
          getBlendOperation(gl, extensions, alpha.operation),
        );
      } else {
        state.disable(gl.BLEND);
      }

      if (colorWrite) {
        state.colorMask(colorWrite[0], colorWrite[1], colorWrite[2], colorWrite[3]);
      } else {
        state.colorMask(true, true, true, true);
      }

      if (depth) {
        state.enable(gl.DEPTH_TEST);
        state.depthFunc(getCompareFunction(gl, depth.compare));
        state.depthMask(depth.write && !depthReadOnly);
      } else {
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);
      }

      if (stencil) {
        state.enable(gl.STENCIL_TEST);
        const { front, back, readMask, writeMask, reference } = stencil;
        state.stencilFuncSeparate(gl.FRONT, getCompareFunction(gl, front.compare), reference, readMask);
        state.stencilFuncSeparate(gl.BACK, getCompareFunction(gl, back.compare), reference, readMask);
        state.stencilOpSeparate(
          gl.FRONT,
          getStencilOperation(gl, front.failOp),
          getStencilOperation(gl, front.depthFailOp),
          getStencilOperation(gl, front.passOp),
        );
        state.stencilOpSeparate(
          gl.BACK,
          getStencilOperation(gl, back.failOp),
          getStencilOperation(gl, back.depthFailOp),
          getStencilOperation(gl, back.passOp),
        );
        state.stencilMask(depthReadOnly ? 0 : writeMask);
      } else {
        state.disable(gl.STENCIL_TEST);
      }

      if (cull === 'front' || cull === 'back') {
        state.enable(gl.CULL_FACE);
        state.cullFace(cull === 'front' ? gl.FRONT : gl.BACK);
      } else {
        state.disable(gl.CULL_FACE);
      }
      state.frontFace(frontFace === 'cw' ? gl.CW : gl.CCW);
    },
    setViewport: (x: number, y: number, width: number, height: number) => {
      state.viewport(x, y, width, height);
    },
    setScissor: (x: number, y: number, width: number, height: number) => {
      state.enable(gl.SCISSOR_TEST);
      state.scissor(x, y, width, height);
    },
    disableScissor: () => {
      state.disable(gl.SCISSOR_TEST);
    },
    clear: (values: ClearValues) => {
      checkRenderPass(state, 'clear');
      const { color, depth, stencil } = values;
      if (depthReadOnly && (depth !== undefined || stencil !== undefined)) {
        throw new Error('clearing the depth or stencil of a render pass with a read-only depth attachment');
      }
      let mask = 0;
      if (color) {
        state.colorMask(true, true, true, true);
        if (isWebGL2(gl) && hasIntegerColorAttachment(gl)) {
          intClearValue[0] = color[0];
          intClearValue[1] = color[1];
          intClearValue[2] = color[2];
          intClearValue[3] = color[3];
          gl.clearBufferiv(gl.COLOR, 0, intClearValue);
        } else {
          state.clearColor(color[0], color[1], color[2], color[3]);
          mask |= gl.COLOR_BUFFER_BIT;
        }
      }
      if (depth !== undefined) {
        state.depthMask(true);
        state.clearDepth(depth);
        mask |= gl.DEPTH_BUFFER_BIT;
      }
      if (stencil !== undefined) {
        state.stencilMask(0xff);
        state.clearStencil(stencil);
        mask |= gl.STENCIL_BUFFER_BIT;
      }
      if (mask) gl.clear(mask);

      if (color) {
        const cw = current.colorWrite;
        if (cw) state.colorMask(cw[0], cw[1], cw[2], cw[3]);
      }
      if (depth !== undefined) state.depthMask((current.depth?.write ?? false) && !depthReadOnly);
      if (stencil !== undefined && current.stencil) state.stencilMask(depthReadOnly ? 0 : current.stencil.writeMask);
    },
  };
}
