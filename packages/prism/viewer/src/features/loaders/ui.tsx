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
import { EntryFormParams } from './operations.js';

// Home's "Download by ID" (parity.md §1, Home's loading): the source and the id, and Download through `load.entry` as
// the person. Home's registry actions already list Mol*'s own Download Structure (every source with its options), so
// this stays minimal: the quick path for an id. The status line says what happened (the viewer posts every person's
// run), and a loaded id answers "already loaded" rather than a second copy; the registry action loads it again.

export function loadersUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: 'loaders.section.entry',
        slot: 'home.import',
        title: 'Download by ID',
        icon: 'download',
        order: 0,
        keywords: ['pdb', 'alphafold', 'afdb', 'emdb', 'model archive', 'fetch', 'entry'],
        render: domView(() => {
          const form = new ParamsFormModel(EntryFormParams);
          const busy = observableState(false);
          const download = async () => {
            const { source, id } = form.values.value;
            if (busy.value || !id.trim() || !ctx.api.load.entry) return;
            busy.set(true);
            try {
              await ctx.api.load.entry({ source, id });
            } finally {
              busy.set(false);
            }
          };
          return (
            <div style="display: grid; gap: var(--msp-space-3)">
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
