/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `test:perf`, the quick allocation tier of the Mol* adapter's hot paths (src/sequence/molstar/_test/*.quick.ts, against
// its `_test/perf-baseline.json`): Seq*'s tier (`perfTests`, packages/sequences/vitest.shared.ts and its AGENTS.md,
// Tests), rooted at the adapter's folder. Opt in, never in CI.
import { join } from 'node:path';
import { defineConfig } from 'vitest/config';
import { perfTests } from '../../sequences/vitest.shared.js';

const root = join(import.meta.dirname, 'src/sequence/molstar');
export default defineConfig({ ...perfTests(root), root });
