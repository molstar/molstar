/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { HistoryStep, ReadonlyState } from '@molstar/prism-ui';
import { ApiArea, ApiResult } from './area.js';

// `history` (operations.md §2, sequence.md §7.3): undo the latest step, or an agent's turn. Not history steps
// themselves. No redo yet.

export class HistoryApi extends ApiArea {
  readonly id = 'history';

  /** The undoable steps, newest first. */
  get steps(): ReadonlyState<readonly HistoryStep[]> {
    return this.ctx.managers.history.steps;
  }

  readonly undo = this.operation({
    id: 'history.undo',
    title: 'Undo',
    params: {},
    history: false,
    when: (ctx) => ctx.managers.history.undoBlocked.value ?? true,
    async run(_, ctx) {
      const r = await ctx.managers.history.undo();
      if (!r.ok) return ApiResult.fail('blocked', r.reason);
      return ApiResult.ok(`Undone: ${r.step.label}.`, { view: true });
    },
  });

  readonly undoTurn = this.operation({
    id: 'history.undo-turn',
    title: "Undo the assistant's turn",
    params: {
      turn: PD.Text('', { description: 'The turn, by id' }),
      discardLater: PD.Boolean(false, { description: 'Also undo the steps made after it' }),
    },
    history: false,
    async run({ turn, discardLater }, ctx) {
      if (!turn) return ApiResult.fail('invalid', 'Which turn?', 'Give the turn id.');
      const r = await ctx.managers.history.undoTurn(turn, { discardLater });
      if (!r.ok) {
        return ApiResult.fail(
          r.laterSteps ? 'confirm' : 'blocked',
          r.reason,
          r.laterSteps ? 'Ask first, then call again with discardLater.' : undefined,
        );
      }
      return ApiResult.ok("Undid the assistant's turn.", { view: true });
    },
  });

  readonly list = this.operation({
    id: 'history.list',
    title: 'List the undoable steps',
    params: {},
    read: true,
    async run(_, ctx) {
      const steps = ctx.managers.history.steps.value;
      if (!steps.length) return ApiResult.ok('Nothing to undo.', {}, { value: steps });
      return ApiResult.ok(
        steps.map((s, i) => `${i + 1}. ${s.label}${s.unsafe ? ` (${s.unsafe})` : ''}`).join('\n'),
        {},
        { value: steps },
      );
    },
  });
}
