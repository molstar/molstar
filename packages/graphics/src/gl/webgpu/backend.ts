/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsBackendDescriptor, GraphicsContextAttribs } from '../backend.js';
import { createWebGPUContext, requestWebGPUDevice } from './context.js';

/** `strict`: see `WebGPUContextOptions` */
export function createWebGPUBackend(options: { strict?: boolean } = {}): GraphicsBackendDescriptor {
  return {
    id: 'webgpu',
    isAvailable: () => typeof navigator !== 'undefined' && !!navigator.gpu,
    create: async (canvas: HTMLCanvasElement, attribs: GraphicsContextAttribs) => {
      const device = await requestWebGPUDevice(attribs.powerPreference);
      if (!device) return undefined;
      return createWebGPUContext(device, {
        width: canvas.width,
        height: canvas.height,
        canvas,
        pixelScale: attribs.pixelScale,
        strict: options.strict,
        // a lost device is replaced by one of the same adapter preference
        requestDevice: () => requestWebGPUDevice(attribs.powerPreference),
      });
    },
  };
}

export const WebGPUBackend = createWebGPUBackend();
