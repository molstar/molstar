/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// Seq*'s own vitest configs (RB0; packages/sequences/AGENTS.md, Tests; the lanes: the root AGENTS.md, Tests): its
// heavy lane (`test:heavy`, also in the root's `pnpm test:heavy`), the browser suite (and its heap guards, `HEAP=1`) and
// the bench. The heap guards and the quick allocation tier every package shares (`heapTests`, `perfTests`) are
// `@molstar/test-harness/vitest`'s. The quick lane is the repository's one run (`pnpm test:unit
// packages/sequences/<package>`, the root `vitest.config.ts`), which reads the same rules; `nodeTests` still serves the
// heavy lane and the composer proof. They return plain objects: this folder has no vitest of its own to import from.

import { globSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { BROWSER_DIRS, GUARD_FLAGS, SOURCES } from '@molstar/test-harness/vitest';

const BROWSER_PATH = /(^|\/)_test\/(.*\/)?(canvas|dom|react)\//;
/** A file that mocks or resets modules needs a module registry of its own (`isolate`). */
const MOCKS = /\bvi\.(mock|doMock|unmock|resetModules)\(/;

/**
 * The lanes' tags, as the root config declares them (Vitest's tags are strict). A test is heavy by nature (a wall-clock
 * timing, a fuzz or oracle run of thousands, a million rows, a large fixture, Chrome or a CLI in a process of its own)
 * when it carries `{ tags: ["heavy"] }`, or its file `// @module-tag heavy`; its quick representative (a short fuzz, a
 * few stories, a small scale) stays untagged.
 */
const TAGS = [
  { name: 'medium', description: 'larger fixtures, sweeps: the root `pnpm test:medium`' },
  {
    name: 'heavy',
    description: 'heavy by nature (timing, fuzz, scale, Chrome, CLIs): `test:heavy`, nightly CI',
    timeout: 120_000,
  },
];

/**
 * The Node suites, with no DOM (browser tests run in `browserTests`). `lane`: `"quick"` is the default, the one every
 * push runs: fast, deterministic, no wall-clock timing, no Chrome, no heap guards; `"heavy"` is `test:heavy`, the tests
 * tagged heavy, run by hand, at checkpoints and nightly; `"all"` every lane. Files share their worker's module
 * registry (`isolate: false`: each file imported the whole package again, most of the run's time); the files that mock
 * or reset modules, found by reading them, run isolated in a project of their own. `dir` is the package's folder.
 */
export const nodeTests = (dir: string, lane: 'quick' | 'heavy' | 'all' = 'quick') => {
  const exclude = [...BROWSER_DIRS, '**/node_modules/**'];
  const mocking = globSync('_test/**/*.test.ts', { cwd: dir })
    .filter((f) => !BROWSER_PATH.test(f) && MOCKS.test(readFileSync(join(dir, f), 'utf8')))
    .sort();
  // A project extending the root concatenates its lists with the root's: the root names no files of its own. Vitest
  // queues projects by name: "isolated" first, so its fresh workers start with the run, not after it.
  const shared = {
    extends: true,
    test: { name: 'shared', include: ['_test/**/*.test.ts'], exclude: [...exclude, ...mocking], isolate: false },
  };
  const isolated = { extends: true, test: { name: 'isolated', include: mocking } };
  return {
    ...SOURCES,
    test: {
      passWithNoTests: true,
      ...(process.env.PROF_SETUP ? { setupFiles: [process.env.PROF_SETUP] } : {}),
      tags: TAGS,
      tagsFilter: lane === 'quick' ? ['!medium', '!heavy'] : lane === 'heavy' ? ['heavy'] : undefined,
      projects: mocking.length > 0 ? [shared, isolated] : [shared],
    },
  };
};

/**
 * Chromium on the installed Chrome (no browser download). `HEAP=1` runs the `*.heap.ts` guards instead, with `gc`
 * exposed, the fixed young generation and precise `performance.memory` (_test/canvas/heap-browser.ts).
 * `playwright` is the package's own `@vitest/browser-playwright` export.
 */
export const browserTests = <Provider>(playwright: (o: object) => Provider) => {
  const heap = process.env.HEAP === '1';
  return {
    ...SOURCES,
    test: {
      include: heap ? ['_test/**/{canvas,dom,react}/**/*.heap.ts'] : ['_test/**/{canvas,dom,react}/**/*.test.{ts,tsx}'],
      passWithNoTests: true,
      fileParallelism: false,
      browser: {
        enabled: true,
        headless: true,
        screenshotFailures: false,
        provider: playwright({
          launchOptions: {
            channel: 'chrome',
            args: heap ? [`--js-flags=${GUARD_FLAGS.join(' ')}`, '--enable-precise-memory-info'] : [],
          },
        }),
        instances: [{ browser: 'chromium' as const, viewport: { width: 1200, height: 900 } }],
      },
    },
  };
};

/**
 * `vitest bench`, outside `test`: the perf phases' numbers (`pnpm --filter <package> bench`). A bench is a test of
 * `bench/` that runs the `bench` fixture (Vitest 5; `core/bench/timed.ts`).
 */
export const benchTests = () => ({
  ...SOURCES,
  test: { benchmark: { include: ['bench/**/*.bench.ts'] }, testTimeout: 300_000, hookTimeout: 300_000 },
});
