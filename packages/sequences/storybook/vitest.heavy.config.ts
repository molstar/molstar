/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `test:heavy`: the Node tests tagged heavy (timing, fuzz, scale, Chrome, CLIs), by hand, at checkpoints
// and nightly (the root's `pnpm test:heavy` runs them too).
import { defineConfig } from 'vitest/config';
import { nodeTests } from '../vitest.shared.js';

export default defineConfig(nodeTests(import.meta.dirname, 'heavy'));
