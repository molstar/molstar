// @vitest-environment jsdom
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui/dom/render';
import { flush } from '@molstar/prism-ui/dom/scope';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { it } from 'vitest';
import { checkLinear, scaling } from '../../../../ui/src/_test/scaling.js';
import type { AgentTurn, CommandRun } from '../../model.js';
import { RunList } from '../agent-turn.js';

// A turn's run list while the answer streams grows about linearly in its runs (`test:scaling`, prism-ui's
// _test/scaling.ts): N runs, 50 streamed tokens (the turn's text set each time) and one run's progress. A `find` per
// row and field per token was N² a token.

function turn(runs: number, text: string, progress: number): AgentTurn {
  const list: CommandRun[] = Array.from({ length: runs }, (_, i) => ({
    id: `r${i}`,
    command: 'color',
    title: `Color part ${i}`,
    source: { kind: 'assistant' },
    status: i === runs - 1 ? 'running' : 'done',
    ...(i === runs - 1 ? { progress } : {}),
  })) as CommandRun[];
  return {
    id: 't',
    prompt: 'Color it',
    source: { kind: 'assistant' },
    status: 'working',
    text,
    runs: list,
    checkpoint: false,
    time: new Date(0),
  } as AgentTurn;
}

it('a run list of N runs, an answer streamed', () => {
  const s = scaling(
    {
      setup: (n) => {
        const state = observableState<AgentTurn>(turn(n, '', 0));
        const box = document.createElement('div');
        document.body.append(box);
        const dispose = render(() => <RunList turn={state.asReadonly} />, box);
        flush();
        return { n, state, box, dispose };
      },
      run: ({ n, state }) => {
        let text = '';
        const runs = state.value.runs;
        for (let k = 0; k < 50; k++) {
          text += 'word ';
          // A token: the turn's text, its runs as they were but the last one's progress.
          const last = { ...runs[n - 1], progress: k / 50 };
          state.set({ ...state.value, text, runs: [...runs.slice(0, n - 1), last] });
          flush();
        }
      },
      teardown: ({ box, dispose }) => {
        dispose();
        box.remove();
      },
    },
    200,
  );
  checkLinear('run list: N runs, 50 tokens', s);
});
