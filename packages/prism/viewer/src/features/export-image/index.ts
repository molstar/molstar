/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ApiResult, memberName, type OperationFn } from '../../api/area.js';
import { operationCommand } from '../../commands/operations.js';
import type { ViewerExtension } from '../../spec.js';
import { downloadBlob } from '../../ui/download.js';
import { ExportApi, LoadViewport, type ViewportLoadParams } from './operations.js';

declare module '../../api.js' {
  // The `export` area (operations.md §2): image here; other exports (models, geometry, animation) are extensions'.
  interface ExtensionApis {
    export: ExportApi;
  }
}

declare module '../../api/load.js' {
  // The viewport state's import, in the core's `load` area (present while the feature is in the spec, as by default).
  interface LoadApi {
    readonly viewport?: OperationFn<typeof ViewportLoadParams>;
  }
}

// The `export-image` feature (viewer-spec.md §3.1): the `export` area with `export.image` and the Image section of the
// Export card; the viewport state (`export.viewport`, `load.viewport`) and the camera in the Viewport section (v5's
// screenshot panel, parity.md). The views are a chunk of their own.

export const ExportImageFeature: ViewerExtension = {
  id: 'export-image',
  title: 'Export image',
  builtIn: true,
  api: { export: (ctx) => new ExportApi(ctx) },
  contributes: () => [{ area: 'load', name: memberName(LoadViewport.id), operation: LoadViewport }],
  commands: (api) => {
    const x = api.export;
    if (!x) return [];
    return [
      {
        ...operationCommand(api, x.image, {
          category: 'Export',
          icon: 'camera',
          title: 'Export an image…',
          alias: ['export image', 'screenshot', 'png', 'save image'],
          action: 'Export',
          view: { title: 'Export › Image', card: 'core.overlay.export', section: 'export-image.section.image' },
        }),
        // The operation gives the image's bytes; the command saves the file through the page.
        async run(args) {
          const r = await api.export!.image(args as never);
          if (!r.ok) throw new Error(r.hint ? `${r.message} ${r.hint}` : r.message);
          if (r.value) downloadBlob(r.value.blob, r.value.filename);
          return {
            message: `Saved ${r.value?.filename} (${r.value?.width} × ${r.value?.height}).`,
            ops: `export.image(${JSON.stringify(args)})`,
          };
        },
      },
      operationCommand(api, x.imageCopy, {
        category: 'Export',
        icon: 'copy',
        title: 'Copy an image',
        alias: ['copy image', 'clipboard'],
        action: 'Copy',
        view: { title: 'Export › Image', card: 'core.overlay.export', section: 'export-image.section.image' },
      }),
      {
        ...operationCommand(api, x.viewport, {
          category: 'Export',
          icon: 'camera',
          title: 'Export the viewport state',
          alias: ['viewport state', 'camera state', 'save camera', 'save view'],
          action: 'Export',
        }),
        async run(args) {
          const r = await api.export!.viewport(args as never);
          if (!r.ok) throw new Error(r.hint ? `${r.message} ${r.hint}` : r.message);
          if (r.value) downloadBlob(r.value.blob, r.value.filename);
          return { message: `Saved ${r.value?.filename}.`, ops: `export.viewport(${JSON.stringify(args)})` };
        },
      },
      ...(api.load.viewport
        ? [
            operationCommand(api, api.load.viewport, {
              category: 'Data',
              icon: 'folder-open',
              title: 'Load a viewport state…',
              alias: ['viewport state', 'load camera', 'restore view', 'open viewport'],
              action: 'Load',
            }),
          ]
        : []),
      {
        id: 'export.open',
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
