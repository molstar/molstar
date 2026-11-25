/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

import { checkLinear, scalingAsync } from '@molstar/test-harness/scaling';
import { expect, it } from 'vitest';
import { flush, render } from '../../dom/index.js';
import { HistoryModel, type HistoryOrigin } from '../history-model.js';
import { HistoryView } from '../history-view.js';

// A long agent turn while History is open is linear (history §4): each run publishes the turn's step again, and the
// view lists only its latest runs (`lastRuns`), not every run as rows (R² for R runs: 12.5M rows for 5,000).

const turn: HistoryOrigin = { kind: 'agent', turn: 't1' };

it('an agent turn of R runs with History open', async () => {
  const s = await scalingAsync(
    {
      setup: async (n) => {
        const h = new HistoryModel<number>({ checkpoint: () => 0, restore: async () => {} }, { maxRuns: 1e9 });
        const el = document.createElement('div');
        document.body.append(el);
        const dispose = render(() => <HistoryView history={h} />, el);
        const end = h.beginGroup('Turn', turn);
        return { h, end, n, el, dispose };
      },
      run: async ({ h, n }) => {
        for (let i = 0; i < n; i++) {
          await h.run(`run ${i}`, turn, async () => i);
          flush();
        }
      },
      teardown: ({ h, end, n, el, dispose }) => {
        expect(h.steps.value[0].runCount).toBe(n);
        end();
        dispose();
        el.remove();
        h.dispose();
      },
    },
    1000,
  );
  checkLinear('history view: an agent turn of R runs', s);
});
