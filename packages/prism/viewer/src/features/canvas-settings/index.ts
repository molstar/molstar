/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { OperationFn } from '../../api/area.js';
import { operationCommand, runAsCommand } from '../../commands/operations.js';
import type { ViewerExtension } from '../../spec.js';
import {
  type CameraParams,
  type ClippingParams,
  canvasSettingsOperations,
  type LightingParams,
  type SettingsParams,
  type XrParams,
} from './operations.js';

declare module '../../api/canvas.js' {
  // What the feature adds to the core's `canvas` area (present while it is in the spec, as by default).
  interface CanvasApi {
    readonly settings: OperationFn<typeof SettingsParams>;
    readonly lighting: OperationFn<typeof LightingParams>;
    readonly camera: OperationFn<typeof CameraParams>;
    readonly clipping: OperationFn<typeof ClippingParams>;
    readonly xr: OperationFn<typeof XrParams>;
  }
}

// The `canvas-settings` feature (viewer-spec.md §3.1): the Canvas card and the canvas's lighting,
// illumination, camera, clipping and XR operations, added to the core's `canvas` area. Its views are a chunk of their own.

export const CanvasSettingsFeature: ViewerExtension = {
  id: 'canvas-settings',
  title: 'Canvas settings',
  builtIn: true,
  contributes: () => canvasSettingsOperations(),
  commands: (api) => {
    const canvas = api.canvas;
    return [
      operationCommand(api, canvas.lighting, {
        category: 'Canvas',
        icon: 'lightbulb',
        alias: ['lighting', 'shadows', 'occlusion', 'outline', 'fog'],
      }),
      // Illumination on or off in one go (the toolbar's button): the lighting's switch, the other of what it is.
      {
        id: 'canvas.illumination',
        title: 'Turn illumination on or off',
        category: 'Canvas',
        icon: 'lightbulb',
        alias: ['illumination', 'path tracing', 'gi'],
        when: () => api.describe('canvas.lighting')?.available ?? 'There is no canvas.',
        async run() {
          const r = await api.canvas.get();
          const lit =
            r.ok && (r.value as { lighting?: { illumination?: boolean } } | undefined)?.lighting?.illumination;
          return {
            ...(await runAsCommand(api, 'canvas.lighting', { illumination: { name: lit ? 'off' : 'on' } })),
            undoable: true,
          };
        },
      },
      operationCommand(api, canvas.camera, {
        category: 'Canvas',
        icon: 'camera',
        alias: ['perspective', 'orthographic', 'field of view'],
      }),
      operationCommand(api, canvas.clipping, { category: 'Canvas', icon: 'scan', alias: ['clip', 'slab'] }),
      // Offered only where WebXR is (the operation's `when`).
      ...(canvas.xr
        ? [
            operationCommand(api, canvas.xr, {
              category: 'Canvas',
              icon: 'glasses',
              alias: ['xr', 'vr', 'ar', 'webxr', 'headset'],
              // The command flips it; the operation takes the state wanted (api.md §1).
              toggle: {
                key: 'on',
                current: async () => {
                  const r = await api.canvas.get();
                  return !!(r.ok && (r.value as { xr?: boolean } | undefined)?.xr);
                },
              },
            }),
          ]
        : []),
    ];
  },
  ui: async (ctx) => (await import('./ui.js')).canvasSettingsUI(ctx),
};
