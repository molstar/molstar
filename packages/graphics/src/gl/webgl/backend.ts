/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsBackendDescriptor, GraphicsContextAttribs } from '../backend.js';
import type { GLRenderingContext } from './compat.js';
import { createContext, type WebGLContextProps } from './context.js';

/** Attributes the renderer requires of the drawing buffer */
const RequiredAttribs: WebGLContextAttributes = {
  alpha: true, // the renderer requires an alpha channel
  depth: true, // the renderer requires a depth buffer
  stencil: true, // the renderer requires a stencil buffer
  premultipliedAlpha: true, // the renderer outputs PMA
};

/** A WebGL backend; `options.uniformBlocks`: see `WebGLContextProps` */
export function createWebGLBackend(
  id: 'webgl' | 'webgl2',
  options: Partial<Pick<WebGLContextProps, 'uniformBlocks'>> = {},
): GraphicsBackendDescriptor {
  return {
    id,
    isAvailable: () =>
      typeof (id === 'webgl2' ? globalThis.WebGL2RenderingContext : globalThis.WebGLRenderingContext) !== 'undefined',
    create: async (canvas: HTMLCanvasElement, attribs: GraphicsContextAttribs) => {
      const { powerPreference, failIfMajorPerformanceCaveat, antialias, preserveDrawingBuffer, pixelScale } = attribs;
      const contextAttribs = {
        ...RequiredAttribs,
        powerPreference,
        failIfMajorPerformanceCaveat,
        antialias,
        preserveDrawingBuffer,
      };
      let gl: GLRenderingContext | null = null;
      try {
        gl = canvas.getContext(id, contextAttribs) as GLRenderingContext | null;
        if (!gl && id === 'webgl')
          gl = canvas.getContext('experimental-webgl', contextAttribs) as GLRenderingContext | null;
      } catch {
        gl = null;
      }
      return gl ? createContext(gl, { pixelScale, uniformBlocks: options.uniformBlocks }) : undefined;
    },
  };
}

/** WebGL 2 */
export const WebGL2Backend = createWebGLBackend('webgl2');
/** WebGL 1 with extensions */
export const WebGLBackend = createWebGLBackend('webgl');
