/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `test:heap`, the deep heap guards of the viewer's hot paths (`*.heap.ts` under any `_test/` of src/): the shared tier
// (`heapTests`, @molstar/test-harness/vitest), rooted at src/. Opt in, at checkpoints, never in CI.
import { join } from 'node:path';
import { defineConfig } from 'vitest/config';
import { heapTests } from '@molstar/test-harness/vitest';

const root = join(import.meta.dirname, 'src');
const base = heapTests();
export default defineConfig({ ...base, root, test: { ...base.test, include: ['**/_test/**/*.heap.ts'] } });
