/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Wall-clock bounds, for every package. The fast tier (`test`, the one CI runs) reads no clock: a bound flakes on a
// loaded machine, and CI's is. A test with a bound (a hostile input refused at once, a million rows in a blink) runs
// twice: in the fast tier with its checks and no bound, and in the slow tier (`test:slow`, the `slow` tag) with both.
// A test whose point is the time (a ratio of two sizes, a storm) is tagged slow instead.

import { expect, it } from "vitest";

/** Runs `f` (sync, or a promise awaited) and checks it took under `ms` milliseconds; returns what `f` returned. */
export type Within = <T>(ms: number, f: () => T) => T;

/** No clock: `f`'s result. */
const unbounded: Within = (_, f) => f();

/** `f` under `ms` milliseconds, its result kept. */
const bounded: Within = (ms, f) => {
  const t0 = performance.now();
  const r = f();
  const check = () => expect(performance.now() - t0, `under ${ms} ms`).toBeLessThan(ms);
  if (r instanceof Promise)
    return r.then((v) => {
      check();
      return v;
    }) as typeof r;
  check();
  return r;
};

/**
 * `it(name, body)` twice: in the fast tier with `within` running its function unbounded, and in the slow tier with
 * the bounds checked (named `name` and its bounds). The body bounds only what it times: `within(1000, () => read(b))`.
 */
export function itWithin(name: string, body: (within: Within) => unknown, timeout?: number) {
  it(name, { timeout }, () => body(unbounded));
  it(`${name}: within its time bounds`, { tags: ["slow"], timeout }, () => body(bounded));
}
