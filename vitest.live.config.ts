/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defaultClientConditions, defaultServerConditions } from 'vite';
import { defineConfig } from 'vitest/config';

// Live runs against real models (`pnpm test:live`): opt in, never part of `test:unit` (whose include, `*.test.*`, does
// not match `*.live.ts`). They read their keys from `.env.local` and skip without them; they cost money and reach the
// network.
export default defineConfig({
  resolve: { conditions: ['molstar-src', ...defaultClientConditions] },
  ssr: { resolve: { conditions: ['molstar-src', ...defaultServerConditions] } },
  test: {
    include: ['**/*.live.ts'],
    exclude: ['**/node_modules/**', '**/lib/**', '**/build/**', '.pnpm-store/**', '.claude/**'],
    testTimeout: 180_000,
    fileParallelism: false,
  },
});
