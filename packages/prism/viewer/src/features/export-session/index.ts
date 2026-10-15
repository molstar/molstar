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
import { type SaveParams, type SavedSession, exportSessionOperations } from './operations.js';

declare module '../../api/session.js' {
  // What the feature adds to the core's `session` area (present while it is in the spec, as by default).
  interface SessionApi {
    readonly save?: OperationFn<typeof SaveParams, SavedSession>;
  }
}

// The `export-session` feature (viewer-spec.md §3.1): `session.save` and the "Session & state" section of the Export
// card, which downloads what it gives. Its views are a chunk of their own.

export const ExportSessionFeature: ViewerExtension = {
  id: 'export-session',
  title: 'Export session',
  builtIn: true,
  contributes: () => exportSessionOperations(),
  commands: (api) => {
    const save = api.session.save;
    if (!save) return [];
    // The command saves the file in the page too: the operation itself only makes the bytes.
    return [
      {
        ...operationCommand(api, save, {
          category: 'Data',
          icon: 'save',
          title: 'Save the session…',
          alias: ['save session', 'save state', 'molx', 'molj', 'export session'],
          action: 'Save',
        }),
        async run(args) {
          const { downloadBlob } = await import('./download.js');
          const r = await api.session.save!(args as never);
          if (!r.ok) throw new Error(r.hint ? `${r.message} ${r.hint}` : r.message);
          if (r.value) downloadBlob(r.value.blob, r.value.filename);
          return { message: r.did, ops: `session.save(${JSON.stringify(args)})` };
        },
      },
    ];
  },
  ui: async (ctx) => (await import('./ui.js')).exportSessionUI(ctx),
};
