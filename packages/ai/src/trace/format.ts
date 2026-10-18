/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { TokensView } from './types.js';

/**
 * How the Inspector, the Assistant's card and the History say numbers (trace.md § The Inspector: one usage
 * formatter).
 */

/** A count of tokens: `812`, `3.2k`, `41k`, `1.5M`. */
export function k(n: number): string {
  if (n >= 999_500) return `${Number((n / 1_000_000).toFixed(n >= 9_950_000 ? 0 : 1))}M`;
  return n >= 1000 ? `${(n / 1000).toFixed(n >= 9950 ? 0 : 1)}k` : String(Math.round(n));
}

/** A duration in milliseconds: `340 ms`, `1.8 s`, `42 s`. */
export function ms(t: number): string {
  if (t < 999.5) return `${Math.round(t)} ms`;
  return `${(t / 1000).toFixed(t < 9950 ? 1 : 0)} s`;
}

/** Credits as people read them: `0.4 credits`, `41 credits`. */
export function credits(n: number): string {
  const v = n < 10 ? Math.round(n * 10) / 10 : Math.round(n);
  return `${v} ${v === 1 ? 'credit' : 'credits'}`;
}

/** The share of the input a provider's cache served, in whole percent; undefined when none was. */
export function cachedShare(t: Pick<TokensView, 'input' | 'cached'>): number | undefined {
  return t.cached && t.input > 0 ? Math.round((t.cached / t.input) * 100) : undefined;
}

/**
 * A call's or turn's usage in one line: `3.2k in (85% cached) · 22 out · 0.4 credits`. Tokens written to the cache
 * join the parenthesis (`85% cached, 1.2k written`); a part with nothing to say is left out; counts the harness
 * estimated are marked `~`.
 */
export function formatUsage(t: TokensView): string {
  const share = cachedShare(t);
  const detail = [share !== undefined ? `${share}% cached` : '', t.written ? `${k(t.written)} written` : '']
    .filter(Boolean)
    .join(', ');
  const approx = t.estimated ? '~' : '';
  return [
    `${approx}${k(t.input)} in${detail ? ` (${detail})` : ''}`,
    `${approx}${k(t.output)} out`,
    t.credits !== undefined ? credits(t.credits) : '',
  ]
    .filter(Boolean)
    .join(' · ');
}
