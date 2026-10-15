/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransform } from '@molstar/core/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ApiArea, ApiResult } from './area.js';

// `session` (operations.md §2): `open` (`.molx`, `.molj`) and `clear` now; `save` with `export-session`, the remote
// operations with `remote-states` (plan phase 3). Both are history steps: Undo brings the scene back.

export class SessionApi extends ApiArea {
  readonly id = 'session';

  readonly open = this.operation({
    id: 'session.open',
    title: 'Open a session',
    params: { file: PD.File({ accept: '.molx,.molj,.json', description: 'A session (.molx) or a state (.molj)' }) },
    modules: ['files'],
    approval: 'destructive',
    async run({ file }, ctx) {
      const f = file as unknown as File | null;
      if (!f) return ApiResult.fail('invalid', 'Which session?', 'Pick a .molx or .molj file.');
      let failed: string | undefined;
      const sub = ctx.plugin.events.log.subscribe((e) => {
        if (e.type === 'error') failed = String(e.message);
      });
      try {
        await ctx.plugin.managers.snapshot.open(f);
      } finally {
        sub.unsubscribe();
      }
      if (failed) return ApiResult.fail('failed', `Could not open ${f.name}: ${failed}`);
      return ApiResult.ok(`Opened ${f.name}.`, { view: true, sequence: true });
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
