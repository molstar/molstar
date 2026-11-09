/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// `test:perf`, the graphics' quick allocation tier: the `*.quick.ts` cases under any `_test/` of `src/` (a frame's
// passes, render items' draws, draw-list replay), bytes a step read exactly against `src/_test/perf-baseline.json`
// (`--update` for an intended move). The shared tier (`perfTests`, @molstar/test-harness/vitest). Opt in, by hand;
// never in CI.
import { join } from 'node:path';
import { defineConfig } from 'vitest/config';
import { perfTests } from '@molstar/test-harness/vitest';

const root = join(import.meta.dirname, 'src');
const base = perfTests(root);
export default defineConfig({ ...base, root, test: { ...base.test, include: ['**/_test/**/*.quick.ts'] } });
