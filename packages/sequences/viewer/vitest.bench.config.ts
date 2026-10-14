/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineConfig } from 'vitest/config';
import { benchTests } from '../vitest.shared.js';

export default defineConfig(benchTests());
