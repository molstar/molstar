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
      ['session', 'Session (.molx), with data'],
      ['state', 'State (.molj)'],
    ],
    {
      description:
        'A session (.molx) keeps the loaded data in the file; a state (.molj) only says how to build the scene again',
    },
  ),
  name: PD.Text('', { description: 'The file name, without its extension (empty: from the loaded entries)' }),
  // journal.md §5: opt-in, as prompts and replies are in it.
  record: PD.Boolean(false, {
    label: 'Record',
    description:
      "What was done in this session, the assistant's prompts and replies included (off by default); it opens in History, read-only",
  }),
  // v5's State panel's Save Options (`PluginState.SnapshotParams`): what of the scene goes into the file.
  options: PD.Group(
    {
      data: PD.Boolean(true, { label: 'Data' }),
      camera: PD.Boolean(true, { label: 'Camera' }),
      canvas3d: PD.Boolean(true, { label: 'Canvas settings' }),
      componentManager: PD.Boolean(true, { label: 'Component options' }),
      structureSelection: PD.Boolean(false, { label: 'Selection' }),
      animation: PD.Boolean(true, { label: 'Animation' }),
      behavior: PD.Boolean(false, { label: 'Behaviors' }),
      interactivity: PD.Boolean(true, { label: 'Interactivity' }),
    },
    { label: 'Save options', description: 'What of the scene goes into the file', isExpanded: false },
  ),
};

/** The record's file in a session's zip (journal.md §5). */
export const RecordFile = 'record.jsonl';

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
      async run({ kind, name, options, record }, ctx): Promise<ApiResult<SavedSession>> {
        const plugin = ctx.plugin;
        try {
          // The record beside the state (journal.md §5): `record.jsonl` in a session's zip, a `record` field (its
          // entries) in a state.
          const lines = record ? await ctx.journal.lines() : '';
          const extra = !lines
            ? undefined
            : kind === 'session'
              ? { files: { [RecordFile]: new TextEncoder().encode(lines) } }
              : { fields: { record: lines.split('\n').map((l) => JSON.parse(l) as unknown) } };
          // Without a canvas there is no camera to save (headless).
          // The person's Snapshots list is theirs: the scene goes into the file only (keepList).
          const blob = await plugin.managers.snapshot.serialize({
            keepList: true,
            type: kind === 'session' ? 'molx' : 'molj',
            params: { ...options, camera: options.camera && !!plugin.canvas3d },
            ...(extra ? { extra } : {}),
          });
          const filename = fileName(name, kind, ctx);
          const size =
            blob.size < 1024 * 1024
              ? `${Math.max(1, Math.round(blob.size / 1024))} KB`
              : `${(blob.size / 1048576).toFixed(1)} MB`;
          const withRecord = lines ? ', with its record' : '';
          return ApiResult.ok(`Saved ${filename} (${size}${withRecord}).`, {}, { value: { blob, filename } });
        } catch (e) {
          return ApiResult.fail('failed', `Could not save the ${kind}: ${e instanceof Error ? e.message : String(e)}`);
        }
      },
    },
  },
];
