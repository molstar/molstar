/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

import { checkLinear, scalingAsync } from '@molstar/test-harness/scaling';
import { expect, it } from 'vitest';
import { HistoryModel, type HistoryOrigin } from '../history-model.js';

// A long agent turn grows its one entry linearly (history §2): each run joins the step's runs without a copy of the runs
// before it (a copy each: R² for R runs).

const turn: HistoryOrigin = { kind: 'agent', turn: 't1' };

it('an agent turn of R runs joining one entry', async () => {
  const s = await scalingAsync(
    {
      setup: async (n) => {
        const h = new HistoryModel<number>({ checkpoint: () => 0, restore: async () => {} });
        const end = h.beginGroup('Turn', turn);
        return { h, end, n };
      },
      run: async ({ h, n }) => {
        for (let i = 0; i < n; i++) {
          await h.run(`run ${i}`, turn, async () => i);
          if (h.steps.value[0].id !== 1) throw new Error('one entry');
        }
      },
      teardown: ({ h, end, n }) => {
        expect(h.steps.value[0].runs.length).toBe(n);
        end();
        h.dispose();
      },
    },
    4000,
  );
  checkLinear('history: an agent turn of R runs', s);
});
