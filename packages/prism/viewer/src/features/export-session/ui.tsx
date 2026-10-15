/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ActionFooter, domView, ParamsForm } from '@molstar/prism-ui';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { downloadBlob } from '../../ui/download.js';
import { SaveParams } from './operations.js';

// The "Session & state" section of the Export card: what to save and under which name, and Download, which asks
// `session.save` for the file and saves it through the page.

export function exportSessionUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: 'export-session.section.session',
        slot: 'export',
        title: 'Session & state',
        icon: 'save',
        order: 1,
        keywords: ['save', 'molx', 'molj', 'state', 'download'],
        render: domView(() => {
          const form = new ParamsFormModel(SaveParams);
          const busy = observableState(false);
          const download = async () => {
            if (busy.value || !ctx.api.session.save) return;
            busy.set(true);
            try {
              const r = await ctx.api.session.save(form.values.value);
              if (!r.ok || !r.value) return;
              downloadBlob(r.value.blob, r.value.filename);
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
                busy={busy.asReadonly}
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
