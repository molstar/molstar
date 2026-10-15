/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type ApiContext, ApiResult, type OperationDefinition } from '../../api/area.js';

// The `export-session` feature's operation (operations.md §2, the `session` row: "save arrives with export-session"):
// `session.save` serializes the plugin's state with its snapshot manager and hands back the file's bytes; saving it
// where the person can reach it (a download) is the page's, in the view. It changes nothing in the scene: no history
// step, nothing to undo.

export type SessionKind = 'session' | 'state';

export const SaveParams = {
  kind: PD.Select<SessionKind>(
    'session',
    [
      ['session', 'Session (.molx): the state and the data it uses'],
      ['state', 'State (.molj): the state only'],
    ],
    { description: 'A session keeps the loaded data in the file; a state only says how to build the scene again' },
  ),
  name: PD.Text('', { description: 'The file name, without its extension (empty: from the loaded entries)' }),
};

/** What `session.save` gives. */
export interface SavedSession {
  blob: Blob;
  filename: string;
}

/** A file name from the person's words: no path, no extension of ours twice. */
function fileName(name: string, kind: SessionKind, ctx: ApiContext): string {
  const ext = kind === 'session' ? '.molx' : '.molj';
  const entered = name
    .trim()
    .replace(/[\\/:*?"<>|]+/g, '-')
    .replace(/\.(molx|molj)$/i, '');
  const loaded = ctx.sources.structures.value.map((s) => s.model || s.label).filter(Boolean);
  const base = entered || [...new Set(loaded)].join('-') || 'molstar-session';
  return base + ext;
}

export const exportSessionOperations = (): {
  area: string;
  name: string;
  operation: OperationDefinition<any, any>;
}[] => [
  {
    area: 'session',
    name: 'save',
    operation: {
      id: 'session.save',
      title: 'Save the session',
      params: SaveParams,
      history: false,
      undoable: false,
      when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
      async run({ kind, name }, ctx): Promise<ApiResult<SavedSession>> {
        const plugin = ctx.plugin;
        try {
          // Without a canvas there is no camera to save (headless).
          const blob = await plugin.managers.snapshot.serialize({
            type: kind === 'session' ? 'molx' : 'molj',
            params: { camera: !!plugin.canvas3d },
          });
          const filename = fileName(name, kind, ctx);
          const size =
            blob.size < 1024 * 1024
              ? `${Math.max(1, Math.round(blob.size / 1024))} KB`
              : `${(blob.size / 1048576).toFixed(1)} MB`;
          return ApiResult.ok(`Saved ${filename} (${size}).`, {}, { value: { blob, filename } });
        } catch (e) {
          return ApiResult.fail('failed', `Could not save the ${kind}: ${e instanceof Error ? e.message : String(e)}`);
        }
      },
    },
  },
];
