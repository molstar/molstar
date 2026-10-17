/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `test:heap`, prism-ui's heap guards (packages/prism/AGENTS.md, Tests): the `*.heap.ts` guards under any `_test/` of
// `src/`, heap growth over n and 4n steps with a positive control each. Seq*'s tier (`heapTests`,
// packages/sequences/vitest.shared.ts). Opt in, at checkpoints; never in CI.
import { join } from 'node:path';
import { defineConfig } from 'vitest/config';
import { heapTests } from '../../sequences/vitest.shared.js';

const root = join(import.meta.dirname, 'src');
const base = heapTests();
export default defineConfig({ ...base, root, test: { ...base.test, include: ['**/_test/**/*.heap.{ts,tsx}'] } });
