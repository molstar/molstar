/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Usage } from '@molstar/ai-harness';

/** A model's list price, in dollars per million tokens. */
export interface ModelPrice {
  input: number;
  cachedInput: number;
  /** Cache writes, where the provider bills them (Claude); else as input. */
  cacheWrite?: number;
  output: number;
}

/** The app's prices, by model id: the one table both money (metering) and credits are read from. */
export type PriceTable = Readonly<Record<string, ModelPrice>>;

type Tokens = Pick<Usage, 'uncachedInput' | 'cacheRead' | 'cacheWrite' | 'output'>;

/** What tokens cost in dollars at a model's list price; undefined for a model without one. */
export function dollarsOf(table: PriceTable, model: string, usage: Tokens): number | undefined {
  const p = table[model];
  if (p === undefined) return undefined;
  return (
    (usage.uncachedInput * p.input +
      usage.cacheRead * p.cachedInput +
      usage.cacheWrite * (p.cacheWrite ?? p.input) +
      usage.output * p.output) /
    1e6
  );
}

/**
 * What model calls cost, in credits: a relative unit, not money, read from the same prices, so a turn on one model
 * and one on another compare. A credit is what 1,000 uncached input tokens cost on the `reference` model. Shown,
 * never a limit: budgets are in tokens (tools-and-session.md §5). Undefined for a model without a price (a custom
 * one), or when the reference has none.
 */
export function creditsOf(table: PriceTable, reference: string, model: string, usage: Tokens): number | undefined {
  const credit = table[reference]?.input;
  const dollars = dollarsOf(table, model, usage);
  return credit === undefined || credit <= 0 || dollars === undefined ? undefined : dollars / (credit / 1000);
}
