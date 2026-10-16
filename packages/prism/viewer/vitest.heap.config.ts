/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `test:heap`, the deep heap guards of the Mol* adapter's hot paths (src/sequence/molstar/_test/*.heap.ts): Seq*'s tier
// (`heapTests`, packages/sequences/vitest.shared.ts), rooted at the adapter's folder. Opt in, at checkpoints, never in CI.
import { join } from 'node:path';
import { defineConfig } from 'vitest/config';
import { heapTests } from '../../sequences/vitest.shared.js';

const root = join(import.meta.dirname, 'src/sequence/molstar');
export default defineConfig({ ...heapTests(), root });
