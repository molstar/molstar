/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransform } from '@molstar/core/state';
import type { ApiContext } from './area.js';
import { ApiArea, ApiResult } from './area.js';
import { FileParam, filesOf } from './files.js';

// `session` (operations.md §2): `open` (`.molx`, `.molj`) and `clear` now; `save` with `export-session`, the remote
// operations with `remote-states` (plan phase 3). Both are history steps: Undo brings the scene back.

/** A session file's record (journal.md §5): `record.jsonl` in a session's zip, a `record` field in a state. */
function recordLinesOf(extra: { files?: Record<string, Uint8Array>; fields?: Record<string, unknown> } | undefined) {
  const file = extra?.files?.['record.jsonl'];
  if (file) return new TextDecoder().decode(file);
  const field = extra?.fields?.record;
  return Array.isArray(field) ? field.map((e) => JSON.stringify(e)).join('\n') : '';
}

/** The nodes an opened session could not make: how many, the first one's error, and whether every top one failed. */
function failedNodes(ctx: ApiContext): { count: number; first: string; all: boolean } {
  const state = ctx.plugin.state.data;
  let count = 0;
  let first = '';
  state.cells.forEach((cell) => {
    if (cell.status !== 'error') return;
    count++;
    first ||= cell.errorText ?? 'an error';
  });
  const top = state.tree.children.get(StateTransform.RootRef);
  const all = count > 0 && top.size > 0 && top.toArray().every((ref) => state.cells.get(ref)?.status === 'error');
  return { count, first, all };
}

export class SessionApi extends ApiArea {
  readonly id = 'session';

  readonly open = this.operation({
    id: 'session.open',
    title: 'Open a session',
    params: { file: FileParam({ accept: '.molx,.molj,.json', description: 'A session (.molx) or a state (.molj)' }) },
    modules: ['files'],
    approval: 'destructive',
    async run({ file }, ctx) {
      const f = filesOf(file)[0];
      if (!f) return ApiResult.fail('invalid', 'Which session?', 'Pick a .molx or .molj file.');
      // A lazy extension whose transformers the session uses loads first (a .molj names them; a .molx is zipped).
      if (/\.(molj|json)$/i.test(f.name)) {
        const transformers = [...(await f.text()).matchAll(/"transformer"\s*:\s*"([^"]+)"/g)].map((m) => m[1]);
        await ctx.extensions.ensureFor({ transformers });
      }
      // Only the open itself failing is fatal (it logs why); a node of the session that cannot be made (a server out of
      // reach) is said as a warning beside what did open, and an error logged by other work changes nothing.
      let opened = false;
      let logged: string | undefined;
      const subs = [
        ctx.plugin.managers.snapshot.events.opened.subscribe(() => {
          opened = true;
        }),
        ctx.plugin.events.log.subscribe((e) => {
          if (e.type === 'error') logged = String(e.message);
        }),
      ];
      let extra: Awaited<ReturnType<typeof ctx.plugin.managers.snapshot.open>>;
      try {
        extra = await ctx.plugin.managers.snapshot.open(f);
      } finally {
        for (const s of subs) s.unsubscribe();
      }
      if (!opened) return ApiResult.fail('failed', `Could not open ${f.name}: ${logged ?? 'it is not a session'}`);
      const failed = failedNodes(ctx);
      if (failed.all) {
        return ApiResult.fail(
          'failed',
          `Could not open ${f.name}: nothing in it could be made (${failed.first}).`,
          'A state names its data; data opened from local files is only in a session (.molx).',
        );
      }
      // A record saved with it (journal.md §5): shown in History, read-only, beside the opened scene.
      const lines = recordLinesOf(extra);
      const shown = !!lines && ctx.journal.show(f.name, lines);
      return ApiResult.ok(
        `Opened ${f.name}.${shown ? ' Its record is in History, read-only.' : ''}`,
        { view: true, sequence: true },
        failed.count
          ? {
              warnings: [
                `${failed.count} part${failed.count === 1 ? '' : 's'} of it could not be made: ${failed.first}`,
              ],
            }
          : {},
      );
    },
  });

  readonly clear = this.operation({
    id: 'session.clear',
    title: 'Clear the session',
    params: {},
    approval: 'destructive',
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    async run(_, ctx) {
      const tree = ctx.plugin.state.data.tree;
      const removed = tree.children.get(StateTransform.RootRef).size;
      await ctx.plugin.state.data.build().delete(StateTransform.RootRef).commit();
      // A cleared scene is small: per-step undo is back on (sequence.md §7.4).
      ctx.managers.history.setCheckpointing(true);
      return ApiResult.ok(`Cleared the session (${removed} item${removed === 1 ? '' : 's'}).`, { view: true });
    },
  });
}
