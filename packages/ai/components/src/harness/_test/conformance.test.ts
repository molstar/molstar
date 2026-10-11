/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { runHarnessConformance } from '@molstar/ai-conformance';
import { createMemoryStorage } from '../storage/memory.js';

// The whole suite against the memory backend.
runHarnessConformance({ name: 'memory', storage: () => createMemoryStorage() });
