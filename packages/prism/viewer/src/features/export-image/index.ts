/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ApiResult } from '../../api/area.js';
import { operationCommand } from '../../commands/operations.js';
import type { ViewerExtension } from '../../spec.js';
import { ExportApi } from './operations.js';

declare module '../../api.js' {
  // The `export` area (operations.md §2): image here; other exports (models, geometry, animation) are extensions'.
  interface ExtensionApis {
    export: ExportApi;
  }
}

// The `export-image` feature (viewer-spec.md §3.1): the `export` area with `export.image` and the Image section of the
// Export card. The views are a chunk of their own.

export const ExportImageFeature: ViewerExtension = {
  id: 'export-image',
  title: 'Export image',
  builtIn: true,
  api: { export: (ctx) => new ExportApi(ctx) },
  commands: (api) => {
    const x = api.export;
    if (!x) return [];
    return [
      operationCommand(api, x.image, {
        category: 'Export',
        icon: 'camera',
        title: 'Export an image…',
        alias: ['export image', 'screenshot', 'png', 'save image'],
        action: 'Export',
      }),
      operationCommand(api, x.imageCopy, {
        category: 'Export',
        icon: 'copy',
        title: 'Copy an image',
        alias: ['copy image', 'clipboard'],
        action: 'Copy',
      }),
      {
        id: 'export.show',
        title: 'Show the Image export',
        category: 'Export',
        icon: 'camera',
        alias: ['image export'],
        // Opens the card: a view, not a change.
        async run() {
          const r = await api.view.open({ card: 'core.overlay.export', section: 'export-image.section.image' });
          if (ApiResult.isFailure(r)) throw new Error(r.message);
          return { message: 'Opened the Export card.', ops: 'view.open({ card: "export" })' };
        },
      },
    ];
  },
  ui: async (ctx) => (await import('./ui.js')).exportImageUI(ctx),
};
