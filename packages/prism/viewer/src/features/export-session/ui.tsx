/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, Icon, ParamsForm } from '@molstar/prism-ui';
import { map } from 'rxjs';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { downloadBlob } from '../../ui/download.js';
import { SaveParams, type SessionKind } from './operations.js';

// The "Session & state" section of the Export card, v5's State panel: State (.molj) and Session (.molx) save through
// `session.save` and the page's download, Open reads one (`session.open`); the file's name and the Save options (what of
// the scene goes into the file) folded under them.

export function exportSessionUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: 'export-session.section.session',
        slot: 'export',
        title: 'Session & state',
        icon: 'save',
        order: 1,
        keywords: ['save', 'open', 'molx', 'molj', 'state', 'session', 'download', 'save options'],
        render: domView(() => {
          const { kind: _, ...rest } = SaveParams;
          const form = new ParamsFormModel(rest);
          const busy = observableState<SessionKind | undefined>(undefined);
          const save = async (kind: SessionKind) => {
            if (busy.value || !ctx.api.session.save) return;
            busy.set(kind);
            try {
              const r = await ctx.api.session.save({ kind, ...form.values.value });
              if (!r.ok || !r.value) return;
              downloadBlob(r.value.blob, r.value.filename);
              ctx.managers.notifications.post(r.did, 'success');
            } finally {
              busy.set(undefined);
            }
          };
          const file = (
            <input
              type="file"
              accept=".molx,.molj"
              hidden
              onChange={(e) => {
                const picked = e.currentTarget.files?.[0];
                e.currentTarget.value = '';
                if (picked) void ctx.api.session.open({ file: picked });
              }}
            />
          ) as HTMLInputElement;
          const saving = busy.asReadonly.pipe(map((b) => b !== undefined));
          return (
            <div style="display: grid; gap: var(--msp-space-3); min-width: 0">
              <div style="display: flex; flex-wrap: wrap; gap: var(--msp-space-2)">
                <Button
                  size="sm"
                  icon={<Icon name="download" />}
                  title="Save how to build the scene again (.molj): the data load again from where they came from; not for local files"
                  disabled={saving}
                  onClick={() => void save('state')}
                >
                  State
                </Button>
                <Button
                  size="sm"
                  icon={<Icon name="download" />}
                  title="Save the scene with its data (.molx)"
                  disabled={saving}
                  onClick={() => void save('session')}
                >
                  Session
                </Button>
                <Button
                  size="sm"
                  icon={<Icon name="folder-open" />}
                  title="Open a saved state (.molj) or session (.molx)"
                  onClick={() => file.click()}
                >
                  Open
                </Button>
                {file}
              </div>
              <p style="margin: 0; color: var(--msp-color-text-muted); font-size: var(--msp-font-size-sm)">
                Saved states and sessions might not open in a future version.
              </p>
              <ParamsForm form={form} />
            </div>
          );
        }),
      },
    ],
  };
}
