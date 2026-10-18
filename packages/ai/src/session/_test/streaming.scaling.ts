/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expect, it } from 'vitest';
import { fixture } from './helpers.js';
import {
  checkLinear,
  quadraticControl,
  scaling,
  scalingAsync,
} from '../../../../sequences/core/_test/model/scaling.js';

// An answer streamed (`test:scaling`, `_test/scaling.ts`): a turn whose one reply arrives in N deltas, the host's
// `onText` hearing the whole text so far at each. The session must do O(delta) work for a delta (session.md §3): a
// `trim()` of the growing text flattened its rope per delta, the answer's length squared.

const WORD = 'word '.repeat(20);

it('an answer of N deltas streamed through a turn', async () => {
  const s = await scalingAsync(
    {
      setup: (n) => {
        const deltas = Array.from({ length: n }, () => WORD);
        const { session } = fixture([{ deltas, text: undefined }]);
        let heard = 0;
        return { n, session, heard: () => heard, hear: () => heard++ };
      },
      run: async ({ n, session, hear, heard }) => {
        // The host keeps only the length, as a view that renders the text would: what it does with the text is its own.
        let length = 0;
        const r = await session.turn({ prompt: 'say it' }, { onText: (t) => ((length = t.length), hear()) });
        expect(r.status).toBe('completed');
        expect(length).toBe(n * WORD.length);
        expect(heard()).toBe(n);
        expect(r.text.length).toBe(n * WORD.length);
      },
    },
    1000,
    { reps: 3, alloc: true },
  );
  checkLinear('session: an answer of N deltas', s);
});

it('the quadratic control is caught', async () => {
  checkLinear('control', await scaling(quadraticControl, 2000), false);
});
