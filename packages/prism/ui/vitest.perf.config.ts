/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `test:perf`, prism-ui's quick allocation tier (packages/prism/AGENTS.md, Tests): the `*.quick.ts` cases under any
// `_test/` of `src/`, bytes a step read exactly against `src/_test/perf-baseline.json` (`--update` for an intended
// move). The shared tier (`perfTests`, @molstar/test-harness/vitest). Opt in, by hand; never in CI.
import { join } from 'node:path';
import { defineConfig } from 'vitest/config';
import { perfTests } from '@molstar/test-harness/vitest';

const root = join(import.meta.dirname, 'src');
const base = perfTests(root);
export default defineConfig({ ...base, root, test: { ...base.test, include: ['**/_test/**/*.quick.ts'] } });
