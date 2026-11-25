/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's quick allocation tier, composed as the stories compose it: the `*.quick.ts` cases under any `_test/`
// of `src/` (a warm emphasis frame of a viewer), bytes a step read exactly (`--update` writes `src/_test/perf-baseline.json`). The shared tier (`perfTests`,
// @molstar/test-harness/vitest). Opt in, by hand, never in CI:
// `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.
import { join } from 'node:path';
import { defineConfig } from 'vitest/config';
import { perfTests } from '@molstar/test-harness/vitest';

const root = join(import.meta.dirname, 'src');
const base = perfTests(root);
export default defineConfig({ ...base, root, test: { ...base.test, include: ['**/_test/**/*.quick.ts'] } });
