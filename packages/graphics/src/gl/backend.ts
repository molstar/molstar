/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsBackend } from './capabilities.js';
import type { GraphicsContext } from './context.js';

/** Context creation options, backends ignore what doesn't apply to them */
export interface GraphicsContextAttribs {
  readonly powerPreference: 'default' | 'high-performance' | 'low-power';
  readonly failIfMajorPerformanceCaveat: boolean;
  readonly antialias: boolean;
  readonly preserveDrawingBuffer: boolean;
  readonly pixelScale: number;
}

/**
 * A render backend that consumers pass in (instead of neutral code importing backend modules), so builds include
 * only the backends their consumers list
 */
export interface GraphicsBackendDescriptor {
  readonly id: GraphicsBackend;
  /** Cheap check without creating a context */
  isAvailable(): boolean;
  /** Creates a context on `canvas`, undefined if the context or device can't be created */
  create(canvas: HTMLCanvasElement, attribs: GraphicsContextAttribs): Promise<GraphicsContext | undefined>;
}

/**
 * Context of the first backend of `backends` that is available and can be created on `canvas`. A canvas keeps the
 * first context type it gave out, so a failed backend may leave later ones of another type unable to create one.
 */
export async function createGraphicsContext(
  canvas: HTMLCanvasElement,
  backends: readonly GraphicsBackendDescriptor[],
  attribs: GraphicsContextAttribs,
): Promise<{ context: GraphicsContext; backend: GraphicsBackendDescriptor } | undefined> {
  for (const backend of backends) {
    if (!backend.isAvailable()) continue;
    const context = await backend.create(canvas, attribs);
    if (context) return { context, backend };
  }
  return undefined;
}
