/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type BenchFn, type BenchRunOptions, test } from 'vitest';

/** A benchmark under Vitest 5 (`bench` is a test's fixture there, no longer an export): `bench(name, fn, options)` of Vitest 4. */
export const timed = (name: string, fn: BenchFn, options?: BenchRunOptions) =>
  test(name, async ({ bench }) => {
    await bench(name, fn).run(options);
  });
