/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of hover-card.heap.ts (`test:perf`): a hover version bumped with the pointer on an item whose card is
// open, the guard's step (hover-card-cases.ts), with its positive control.

import { expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { openCard } from './hover-card-cases.js';

it("a hover change on the open card's item", async ({ task }) => {
  const { card, step } = await openCard();
  expect(card.current.value).not.toBe(null);
  quickCase(task, "hover-card: a hover change on the open card's item", step);
});
