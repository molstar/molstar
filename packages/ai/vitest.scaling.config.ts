/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `test:scaling`, the AI package's scaling guards: the `*.scaling.ts` cases under any `_test/` of `src/`, each at n and
// 4n, time and allocation about linear, with a quadratic positive control (`src/session/_test/scaling.ts`). The heap
// tier's process (one thread, `--expose-gc`), as the Prism packages have it. Opt in, at checkpoints; never in CI.
import { join } from 'node:path';
import { defineConfig } from 'vitest/config';
import { heapTests } from '../sequences/vitest.shared.js';

const root = join(import.meta.dirname, 'src');
const base = heapTests();
export default defineConfig({
  ...base,
  root,
  test: { ...base.test, include: ['**/_test/**/*.scaling.ts'], testTimeout: 300_000 },
});
