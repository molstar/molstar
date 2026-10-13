/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defaultClientConditions, defaultServerConditions } from 'vite';
import { defineConfig } from 'vitest/config';

const include = ['**/*.test.{ts,tsx,mjs}'];
const exclude = ['**/node_modules/**', '**/lib/**', '**/build/**', '.pnpm-store/**', '.claude/**', 'smoke/**'];

// The AI packages' tests (the harness conformance suite above all) count and wait for rounds of real work (timers,
// WebCrypto, fake IndexedDB), and lose those races when other test files share the machine: their files run one by
// one, after everything else (as they did when they had a Jest run of their own).
const serial = ['packages/ai/**/*.test.ts'];

// One Vitest run for the whole workspace, ESM throughout. Workspace packages resolve to their sources through the
// `molstar-src` export condition (the same one the build and the type checker use), so a new package needs no mapping
// here. DOM tests opt in per file with `// @vitest-environment jsdom`.
export default defineConfig({
  resolve: { conditions: ['molstar-src', ...defaultClientConditions] },
  ssr: { resolve: { conditions: ['molstar-src', ...defaultServerConditions] } },
  test: {
    // Preact (the adapter tests) ships ESM and CommonJS builds; inlined, all of it comes from one.
    server: { deps: { inline: ['preact'] } },
    // The long cases (2,000-append histories, fresh module graphs) can share the machine with other files.
    testTimeout: 30_000,
    hookTimeout: 30_000,
    environmentOptions: { jsdom: { url: 'http://localhost/' } },
    projects: [
      { extends: true, test: { name: 'unit', include, exclude: [...exclude, ...serial], sequence: { groupOrder: 0 } } },
      {
        extends: true,
        test: { name: 'serial', include: serial, exclude, fileParallelism: false, sequence: { groupOrder: 1 } },
      },
    ],
  },
});
