/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { creditsOf, dollarsOf, type PriceTable } from '../index.js';

// The app's prices (dollars per million tokens); the package has no table of its own.
const table: PriceTable = {
  flash: { input: 0.75, cachedInput: 0.075, output: 3.75 },
  sonnet: { input: 2, cachedInput: 0.2, cacheWrite: 2.5, output: 10 },
};
const usage = { uncachedInput: 10_000, cacheRead: 90_000, cacheWrite: 0, output: 1_000 };

describe('prices (tools-and-session.md §5)', () => {
  it('prices tokens in dollars; a model without a price has none', () => {
    expect(dollarsOf(table, 'flash', usage)).toBeCloseTo((10_000 * 0.75 + 90_000 * 0.075 + 1_000 * 3.75) / 1e6);
    // Cache writes are billed as input where the model has no rate for them, else at theirs.
    expect(dollarsOf(table, 'flash', { ...usage, cacheWrite: 1_000_000 })).toBeCloseTo(0.75 + 0.018, 4);
    expect(dollarsOf(table, 'sonnet', { ...usage, cacheWrite: 1_000_000 })).toBeCloseTo(2.5 + 0.048, 4);
    expect(dollarsOf(table, 'custom', usage)).toBeUndefined();
  });

  it('prices credits relative to the reference model: 1,000 uncached input tokens on it', () => {
    expect(creditsOf(table, 'flash', 'flash', usage)).toBeCloseTo(10 + 9 + 5);
    expect(creditsOf(table, 'flash', 'sonnet', usage)).toBeCloseTo(26.67 + 24 + 13.33, 1);
    expect(creditsOf(table, 'flash', 'a-custom-model', usage)).toBeUndefined();
    expect(creditsOf(table, 'no-reference', 'flash', usage)).toBeUndefined();
  });
});
