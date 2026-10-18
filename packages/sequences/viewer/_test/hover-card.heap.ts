/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover card's hover handler (seqspec:interaction#the-hover-card) at pointer rate: the hover's version bumped with
// the pointer on an item whose card is open (what each move within the item does to the handler) allocates nothing, on
// RB0b's stepped guard; the card and the step are hover-card-cases.ts's, shared with the quick twin
// (hover-card.quick.ts). The positive control: the same bumps allocating an object each (`leaking`), measured alike,
// are caught.

import { expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';
import { N, openCard } from './hover-card-cases.js';

it("a hover change on the open card's item allocates nothing", async () => {
  const { card, step } = await openCard();
  expect(card.current.value).not.toBe(null);
  const r = heapScaling(N, step);
  expect(isFlat(r), describeScaling(r, N)).toBe(true);
});

it('the same bumps allocating an object each are caught (positive control)', async () => {
  const { step } = await openCard();
  expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
});
