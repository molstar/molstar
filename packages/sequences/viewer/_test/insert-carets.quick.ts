/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of insert-carets.heap.ts (`test:perf`): the same hover sweep (insert-carets-cases.ts), a hover and an
// overlay frame of both sections a step, with its positive control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { insertCarets } from './insert-carets-cases.js';

describe('carets over a hover sweep, quick', () => {
  it('a hover and an overlay frame of both sections a move', async ({ task }) => {
    const { frame, share } = await insertCarets();
    quickCase(task, 'insert-carets: hover sweep', frame);
    // Not vacuous: 8 of the 12 points light a caret.
    expect(share()).toBeGreaterThan(0.5);
  });
});
