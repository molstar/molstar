/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { OperationFn } from '../../api/area.js';
import { operationCommand } from '../../commands/operations.js';
import type { ViewerExtension } from '../../spec.js';
import {
  type CameraParams,
  type ClippingParams,
  canvasSettingsOperations,
  type IlluminationParams,
  type LightingParams,
  type SettingsParams,
} from './operations.js';

declare module '../../api/canvas.js' {
  // What the feature adds to the core's `canvas` area (present while it is in the spec, as by default).
  interface CanvasApi {
    readonly settings?: OperationFn<typeof SettingsParams>;
    readonly lighting?: OperationFn<typeof LightingParams>;
    readonly illumination?: OperationFn<typeof IlluminationParams>;
    readonly camera?: OperationFn<typeof CameraParams>;
    readonly clipping?: OperationFn<typeof ClippingParams>;
  }
}

// The `canvas-settings` feature (viewer-spec.md §3.1, plan 3.2's first): the Canvas card and the canvas's lighting,
// illumination, camera and clipping operations, added to the core's `canvas` area. Its views are a chunk of their own.

export const CanvasSettingsFeature: ViewerExtension = {
  id: 'canvas-settings',
  title: 'Canvas settings',
  builtIn: true,
  contributes: () => canvasSettingsOperations(),
  commands: (api) => {
    const canvas = api.canvas;
    if (!canvas.lighting || !canvas.illumination || !canvas.camera || !canvas.clipping) return [];
    return [
      operationCommand(api, canvas.lighting, {
        category: 'Canvas',
        icon: 'lightbulb',
        alias: ['lighting', 'shadows', 'occlusion', 'outline', 'fog'],
      }),
      operationCommand(api, canvas.illumination, {
        category: 'Canvas',
        icon: 'lightbulb',
        alias: ['illumination', 'path tracing', 'gi'],
      }),
      operationCommand(api, canvas.camera, {
        category: 'Canvas',
        icon: 'camera',
        alias: ['perspective', 'orthographic', 'field of view'],
      }),
      operationCommand(api, canvas.clipping, { category: 'Canvas', icon: 'scan', alias: ['clip', 'slab'] }),
    ];
  },
  ui: async (ctx) => (await import('./ui.js')).canvasSettingsUI(ctx),
};
