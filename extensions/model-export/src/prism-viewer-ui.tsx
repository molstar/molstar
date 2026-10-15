/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ActionFooter, domView, mapState, observableState, ParamsForm } from '@molstar/prism-ui';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import { downloadBlob } from '@molstar/prism-viewer/ui/download';
import { ModelExportParams, SectionId } from './prism-viewer.js';

// The Models section of the Export card (parity.md §6, the mock's look): the format and "only the selection", and
// Download, which asks `model-export.export` for the file (every loaded structure; a ZIP when several) and saves it
// through the page. It acts only through the person's API.

const FormParams = { format: ModelExportParams.format, selection: ModelExportParams.selection };

export function modelExportUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: SectionId,
        slot: 'export',
        title: 'Models',
        icon: 'file',
        order: 12,
        keywords: ['mmcif', 'bcif', 'cif', 'pdb', 'coordinates', 'structure'],
        render: domView(() => {
          const form = new ParamsFormModel(FormParams);
          const busy = observableState(false);
          const count = mapState(ctx.sources.structures, (s) => s.length);
          const download = async () => {
            const op = ctx.api['model-export']?.export;
            if (busy.value || !op) return;
            busy.set(true);
            try {
              const r = await op(form.values.value);
              if (!r.ok) {
                ctx.managers.notifications.post(r.hint ? `${r.message} ${r.hint}` : r.message, 'error');
                return;
              }
              if (r.value) downloadBlob(r.value.blob, r.value.filename);
              ctx.managers.notifications.post(r.did, 'success');
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
                info={mapState(count, (n) => (n > 1 ? `${n} structures, as a ZIP` : ''))}
                busy={busy.asReadonly}
                disabled={mapState(count, (n) => n === 0)}
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
