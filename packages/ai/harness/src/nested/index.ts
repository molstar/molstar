/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '@molstar/ai-agent';
import { HarnessError } from '../contracts/harness.js';
import { MAIN } from '../contracts/records.js';
import type { RunResult, RunSpec } from '../contracts/runner.js';
import { modelStep } from '../loop/model.js';
import { append, boundary, linked, type Run, type Via } from '../loop/run.js';
import { Exit, local, type Staged, ZERO } from '../loop/state.js';
import { planCalls, toolCalls } from '../tools/plan.js';
import { toolsStep } from '../tools/step.js';

/**
 * A nested run inside a call (harness.md §9): the same engine, in memory, staged into the parent's step.
 * It ends `completed` with its last message, or `failed` when a decision answers `fail`; a stop or a
 * park ends the whole step.
 */
export async function nest(parent: Run, spec: RunSpec, via: Via): Promise<RunResult> {
  const depth = parent.depth + 1;
  if (depth > parent.core.maxDepth)
    return { status: 'failed', usage: ZERO, failure: { reason: 'depth', limit: parent.core.maxDepth } };
  if (spec.pointer === MAIN) throw new HarnessError('pointer-busy', 'a nested run cannot write main');
  const { store } = parent.core;
  let context: Staged;
  if (spec.pointer !== undefined) context = via.journal.open(spec.pointer);
  else if (spec.from !== undefined) {
    const entries = store.heads.get(spec.from);
    if (entries === undefined) throw new HarnessError('invalid', `no head ${spec.from}`);
    context = local(entries);
  } else context = local();
  const abort = linked(via.signal);
  const lastCompaction = spec.pointer !== undefined ? store.lastCompaction(spec.pointer) : undefined;
  const run: Run = {
    core: parent.core,
    ...(parent.turnId !== undefined ? { turnId: parent.turnId } : {}),
    step: parent.step,
    depth,
    path: [...parent.path, via.callId],
    parentCallId: via.callId,
    parentKey: via.key,
    model: spec.model,
    tools: spec.tools,
    definitions: spec.tools.definitions,
    budget: { limits: spec.limits, consumed: ZERO },
    drives: spec.drive !== undefined ? [spec.drive, ...parent.drives] : parent.drives,
    context,
    collector: parent.collector,
    state: parent.state,
    abort,
    ...(parent.configuration !== undefined ? { configuration: parent.configuration } : {}),
    seen: new Map(),
    ...(lastCompaction !== undefined ? { lastCompaction } : {}),
  };
  try {
    append(run, context, spec.input);
    for (;;) {
      const message = await modelStep(run);
      if (toolCalls(message).length === 0) {
        append(run, context, [message]);
        await boundary(run, true);
        return { status: 'completed', output: message, usage: run.budget.consumed };
      }
      const { round, ids } = await toolsStep(run, message, planCalls(run, message, run.step));
      append(run, context, round, ids);
      await boundary(run, false);
    }
  } catch (error) {
    if (error instanceof Exit && error.kind === 'fail' && !run.state.exits.includes(error))
      return { status: 'failed', usage: run.budget.consumed, failure: (error.reason ?? null) as Json };
    throw error;
  }
}
