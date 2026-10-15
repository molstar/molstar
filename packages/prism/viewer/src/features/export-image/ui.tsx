/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView, onCleanup } from '@molstar/prism-ui';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { ImageExport } from '../../views/export/image-export.js';
import { createImageExportModel } from './model.js';

// The Image section of the Export card (parity.md §4's export row): the existing view over a model built on
// `export.preview`, `export.image` and `export.image-copy`; the model is made when the section is shown and released
// with it.

export function exportImageUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: 'export-image.section.image',
        slot: 'export',
        title: 'Image',
        icon: 'camera',
        order: 0,
        keywords: ['screenshot', 'png', 'picture', 'copy'],
        render: domView(() => {
          const model = createImageExportModel({
            api: ctx.api,
            sources: ctx.sources,
            notify: (text, kind) => ctx.managers.notifications.post(text, kind),
          });
          onCleanup(() => model.dispose());
          return <ImageExport model={model} />;
        }),
      },
    ],
  };
}
