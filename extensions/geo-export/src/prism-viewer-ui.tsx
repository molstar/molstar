/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ActionFooter, domView, observableState, ParamsForm } from '@molstar/prism-ui';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import { downloadBlob } from '@molstar/prism-viewer/ui/download';
import { GeometryExportParams, SectionId } from './prism-viewer.js';

// The Geometry section of the Export card (parity.md §6, the mock's look): the format, and Download, which asks
// `geo-export.export` for the file and saves it through the page. It acts only through the person's API.

const FormParams = { format: GeometryExportParams.format };

export function geoExportUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: SectionId,
        slot: 'export',
        title: 'Geometry',
        icon: 'box',
        order: 11,
        keywords: ['glb', 'gltf', 'stl', 'obj', 'usdz', 'mesh', '3d print', 'ar'],
        render: domView(() => {
          const form = new ParamsFormModel(FormParams);
          const busy = observableState(false);
          const download = async () => {
            const op = ctx.api['geo-export']?.export;
            if (busy.value || !op) return;
            busy.set(true);
            try {
              // The run path says how it went (a warning on failure, the result when it worked).
              const r = await op(form.values.value);
              if (r.ok && r.value) downloadBlob(r.value.blob, r.value.filename);
            } finally {
              busy.set(false);
            }
          };
          return (
            <div style="display: grid; gap: var(--msp-space-3); min-width: 0">
              <ParamsForm form={form} />
              <ActionFooter
                icon="download"
                label="Download"
                busy={busy.asReadonly}
                disabled={ctx.sources.empty}
                onReset={() => form.reset()}
                onApply={() => void download()}
              />
            </div>
          );
        }),
      },
    ],
  };
}
