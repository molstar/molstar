/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `test:perf`, the quick allocation tier (`perfTests`, ../vitest.shared.ts).
import { defineConfig } from 'vitest/config';
import { perfTests } from '../vitest.shared.js';

export default defineConfig(perfTests(import.meta.dirname));
