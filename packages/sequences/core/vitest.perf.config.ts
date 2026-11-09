/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `test:perf`, the quick allocation tier (`perfTests`, @molstar/test-harness/vitest).
import { defineConfig } from 'vitest/config';
import { perfTests } from '@molstar/test-harness/vitest';

export default defineConfig(perfTests(import.meta.dirname));
