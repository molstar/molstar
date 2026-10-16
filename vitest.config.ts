/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { globSync, readFileSync } from 'node:fs';
import { defaultClientConditions, defaultServerConditions } from 'vite';
import { defineConfig } from 'vitest/config';

const include = ['**/*.test.{ts,tsx,mjs}'];
const exclude = ['**/node_modules/**', '**/lib/**', '**/build/**', '.pnpm-store/**', '.claude/**', 'smoke/**'];

// Seq* (`packages/sequences`, .v6/seq/plan.md §2.3): its fast tier runs here, in two projects of its own. Files share a
// worker's module registry (`isolate: false`: each file importing the whole package again was most of the run's time);
// the files that mock or reset modules, found by reading them, run isolated. Never matched: the browser folders, the
// heap guards (`*.heap.ts`, not `*.test.*`) and the benches. The opt-in tiers (`test:slow`, `test:browser`,
// `test:heap`, `bench`) are the packages' own scripts (`packages/sequences/vitest.shared.ts`).
const seq = 'packages/sequences';
const seqBrowser = [`${seq}/*/_test/**/{canvas,dom,react}/**`, `${seq}/*/bench/**`, `${seq}/**/*.heap.ts`];
const seqMocks = /\bvi\.(mock|doMock|unmock|resetModules)\(/;
const seqTests = globSync(`${seq}/*/_test/**/*.test.ts`)
  .map((f) => f.replaceAll('\\', '/'))
  .filter((f) => !/\/_test\/(.*\/)?(canvas|dom|react)\//.test(f))
  .sort();
const seqMocking = seqTests.filter((f) => seqMocks.test(readFileSync(f, 'utf8')));

// The slow tier's tag (Seq*: `{ tags: ['slow'] }` on a test, `// @module-tag slow` on a file). Vitest's tags are strict,
// so it is declared; `tagsFilter` is CLI-only in Vitest 5, so the tag itself skips its tests here. Only the Seq*
// projects declare it: Mol*'s own projects are as before. The packages' `test:slow` configs run it (§2.3).
const slow = {
  name: 'slow',
  description: 'slow by nature (timing, fuzz, scale, Chrome, CLIs): the package `test:slow`',
  skip: true,
};

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
      {
        extends: true,
        test: { name: 'unit', include, exclude: [...exclude, `${seq}/**`], sequence: { groupOrder: 0 } },
      },
      // Seq*: listed by file (a glob would have to be inverted to leave the mocking files out of this project).
      {
        extends: true,
        test: {
          name: 'sequences',
          include: seqTests.filter((f) => !seqMocking.includes(f)),
          exclude: [...exclude, ...seqBrowser],
          tags: [slow],
          isolate: false,
          sequence: { groupOrder: 0 },
        },
      },
      {
        extends: true,
        test: {
          name: 'sequences-isolated',
          include: seqMocking,
          tags: [slow],
          exclude: [...exclude, ...seqBrowser],
          sequence: { groupOrder: 0 },
        },
      },
    ],
  },
});
