/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineConfig } from 'vitest/config';
import { nodeTests } from '../vitest.shared.js';

// Three workers: each imports Mol*, most of a worker's start; more workers made the run slower (2026-10-14).
const base = nodeTests(import.meta.dirname);
export default defineConfig({ test: { ...base.test, maxWorkers: 3 } });
