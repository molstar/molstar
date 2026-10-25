/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsBackendDescriptor, GraphicsContextAttribs } from '../backend.js';
import { createWebGPUContext, requestWebGPUDevice, type WebGPUFeatures, type WebGPULimits } from './context.js';

/**
 * `strict`: see `WebGPUContextOptions`; `limits`: the adapter's maxima (default) or the core defaults of the spec;
 * `features`: the optional features the adapter has (default) or none; `requiredFeatures`: adapters without them create
 * no context, so that the next backend of a list is used
 */
export function createWebGPUBackend(
  options: {
    strict?: boolean;
    /** the binding limit of the theme layers' pages, below the device's (`WebGPUContextOptions.maxDataBindingSize`) */
    maxDataBindingSize?: number;
    limits?: WebGPULimits;
    features?: WebGPUFeatures;
    requiredFeatures?: readonly GPUFeatureName[];
  } = {},
): GraphicsBackendDescriptor {
  const request = (powerPreference: GraphicsContextAttribs['powerPreference']) =>
    requestWebGPUDevice(powerPreference, undefined, options.limits, options.features, options.requiredFeatures);
  return {
    id: 'webgpu',
    isAvailable: () => typeof navigator !== 'undefined' && !!navigator.gpu,
    create: async (canvas: HTMLCanvasElement, attribs: GraphicsContextAttribs) => {
      const device = await request(attribs.powerPreference);
      if (!device) return undefined;
      return createWebGPUContext(device, {
        width: canvas.width,
        height: canvas.height,
        canvas,
        pixelScale: attribs.pixelScale,
        strict: options.strict,
        maxDataBindingSize: options.maxDataBindingSize,
        // a lost device is replaced by one of the same adapter preference
        requestDevice: () => request(attribs.powerPreference),
        ownsDevice: true,
      });
    },
  };
}

export const WebGPUBackend = createWebGPUBackend();

/**
 * WebGPU for default backend lists, where it competes with WebGL2: skipped (the next backend of the list is used) where
 * the adapter cannot blend into float32 targets (`float32-blendable`, e.g. older Firefox), as WebGL2 then has more
 * features (the transparency modes wboit and dpoit and illumination, with `EXT_float_blend`); `WebGPUBackend` when
 * WebGPU is asked for
 */
export const AutoWebGPUBackend = createWebGPUBackend({ requiredFeatures: ['float32-blendable'] });
