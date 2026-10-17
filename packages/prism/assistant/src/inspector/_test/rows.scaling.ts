/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ToolStepView, TraceData, TurnView } from '@molstar/ai/trace';
import { observableState } from '@molstar/prism-ui';
import type { Observable } from 'rxjs';
import { expect, it } from 'vitest';
import { checkLinear, quadraticControl, scaling } from '../../../../ui/src/_test/scaling.js';
import { MockTools, mockLimits } from '../../testing/mock-inspector.js';
import { InspectorModel, type TraceRow } from '../model.js';

/** The latest value of a stream that replays at once (the model, with no settling). */
function latest(rows: Observable<readonly TraceRow[]>): readonly TraceRow[] {
  let value: readonly TraceRow[] = [];
  rows.subscribe((r) => (value = r)).unsubscribe();
  return value;
}

// The Inspector's rows as the running turn grows (`test:scaling`, prism-ui's _test/scaling.ts; trace.md, Incremental
// reads): N finished turns of 6 steps each and the running one, which a re-read gives with one more step (the turns that
// ended as the same objects). Rows are built for the turn that changed; the finished turns' rows are the ones kept, so
// the work is the list put together again (about the factor, 4×), not a row made again for every step (which would also
// rebuild every view on them: the allocation would show it).

const step = (id: string): ToolStepView => ({
  kind: 'tool',
  id,
  callId: id,
  toolId: 'echo',
  name: 'echo',
  args: { id },
  status: 'succeeded',
  attempts: [{ ordinal: 1, start: 0, end: 5, status: 'succeeded' }],
  result: 'ok',
});

const turn = (id: string, steps: number, status: TurnView['status']): TurnView => ({
  id,
  prompt: `Turn ${id}`,
  status,
  start: 0,
  tokens: { input: 1, output: 1 },
  steps: Array.from({ length: steps }, (_, i) => step(`${id}-s${i}`)),
  milestones: [],
});

const trace = (turns: TurnView[]): TraceData =>
  ({
    harnessId: 'h',
    execution: { seq: 1, headId: 'h' },
    budget: { consumed: 0, ceiling: 0 },
    preamble: { guidance: '', guidanceTokens: 0, tools: [], toolsTokens: 0 },
    loaded: [],
    declared: [],
    turns,
  }) satisfies TraceData;

it('the rows after the running turn grew by a step, N turns before it', () => {
  const s = scaling(
    {
      setup: (n) => {
        const done = Array.from({ length: n }, (_, i) => turn(`t${i}`, 6, 'done'));
        const running = turn('run', 3, 'working');
        const first = trace([...done, running]);
        const state = observableState<TraceData | undefined>(first);
        const model = new InspectorModel({
          trace: state.asReadonly,
          tools: observableState(MockTools).asReadonly,
          limits: observableState(mockLimits(first)).asReadonly,
          settle: 0,
        });
        const before = latest(model.rows);
        const grown = { ...running, steps: [...running.steps, step('run-new')] };
        return { model, state, done, grown, before };
      },
      run: ({ model, state, done, grown, before }) => {
        state.set(trace([...done, grown]));
        const rows = latest(model.rows);
        expect(rows.length).toBe(before.length + 1);
      },
    },
    250,
    { reps: 7 },
  );
  checkLinear('inspector rows: a step added to the running turn, N turns', s);
});

it('the quadratic control is caught', () => {
  checkLinear('control', scaling(quadraticControl, 2000), false);
});
