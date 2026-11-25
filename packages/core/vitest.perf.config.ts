/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Core's quick allocation tier: the `*.quick.ts` cases under any `_test/` of `src/`, bytes a step read exactly (the
// shared tier, `perfTests` in @molstar/test-harness/vitest). The cases state their bounds themselves (no baseline yet).
// Opt in, by hand, never in CI: `npx vitest run -c packages/core/vitest.perf.config.ts`.
import { join } from 'node:path';
import { defineConfig } from 'vitest/config';
import { perfTests } from '@molstar/test-harness/vitest';

const root = join(import.meta.dirname, 'src');
const base = perfTests(root);
// biome-ignore lint/style/noDefaultExport: Vitest reads a config's default export (biome.json lists the other packages' configs)
export default defineConfig({ ...base, root, test: { ...base.test, include: ['**/_test/**/*.quick.ts'] } });
