/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expect, it } from 'vitest';
import { checkLinear, scalingAsync } from '../../_test/scaling.js';
import { HistoryManager, type HistoryOrigin } from '../history.js';

// A long agent turn grows its one entry linearly (`test:scaling`, ui/src/_test/scaling.ts): each run joins the step's
// runs without a copy of the runs before it (a copy each: R² for R runs).

const turn: HistoryOrigin = { kind: 'agent', turn: 't1' };

it('an agent turn of R runs joining one entry', async () => {
  const s = await scalingAsync(
    {
      setup: async (n) => {
        const h = new HistoryManager<number>({ checkpoint: () => 0, restore: async () => {} });
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
