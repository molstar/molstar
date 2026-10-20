/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterAll, vi } from 'vitest';

// Files share their worker's module registry (vitest.config.ts): what a file changed of it is put back when the file
// ends, so the next file in the worker starts as it would alone. Mocked modules are not undone here: a file that mocks
// or resets modules runs in a registry of its own.
afterAll(() => {
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
  vi.unstubAllEnvs();
  vi.useRealTimers();
});
