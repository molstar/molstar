/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The vitest configs the opt-in tiers share (RB0; .v6/seq/plan.md §2.3). A package's own configs call these: the slow
// Node tier (`test:slow`), the heap guards, the browser suite (and its heap guards, `HEAP=1`) and the bench. The fast
// tier (`test`) is the repository's one run (`pnpm test:unit packages/sequences/<package>`, the root `vitest.config.ts`),
// which reads the same rules; `nodeTests` still serves the slow tier and the composer proof. They return plain objects:
// this folder has no vitest of its own to import from.

import { existsSync, globSync, readFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { defaultClientConditions, defaultServerConditions } from 'vite';
import { PerfReporter } from './perf-reporter.js';

const GUARD_FLAGS = ['--expose-gc', '--min-semi-space-size=64', '--max-semi-space-size=64'];
const BROWSER_DIRS = ['_test/**/canvas/**', '_test/**/dom/**', '_test/**/react/**'];
const BROWSER_PATH = /(^|\/)_test\/(.*\/)?(canvas|dom|react)\//;

/** Workspace packages resolve to their sources through the `molstar-src` export condition, as in the root config. */
const SOURCES = {
  resolve: { conditions: ['molstar-src', ...defaultClientConditions] },
  ssr: { resolve: { conditions: ['molstar-src', ...defaultServerConditions] } },
};
/** A file that mocks or resets modules needs a module registry of its own (`isolate`). */
const MOCKS = /\bvi\.(mock|doMock|unmock|resetModules)\(/;

/**
 * The slow tier's tag. A test is slow by nature (a wall-clock timing, a fuzz or oracle run of thousands, a million
 * rows, a large fixture, Chrome or a CLI in a process of its own) when it carries `{ tags: ["slow"] }`, or its file
 * `// @module-tag slow`; its quick representative (a short fuzz, a few stories, a small scale) stays untagged.
 */
const SLOW_TAG = {
  name: 'slow',
  description: 'slow by nature (timing, fuzz, scale, Chrome, CLIs): `test:slow`, never CI',
  timeout: 120_000,
};

/**
 * The Node suites, with no DOM (browser tests run in `browserTests`). `tier`: `"fast"` is `test`, the default and the
 * one CI runs: fast, deterministic, no wall-clock timing, no Chrome, no heap guards; `"slow"` is `test:slow`, the tests
 * tagged slow, run by hand or at checkpoints; `"all"` both. Files share their worker's module registry (`isolate:
 * false`: each file imported the whole package again, most of the run's time); the files that mock or reset modules,
 * found by reading them, run isolated in a project of their own. `dir` is the package's folder.
 */
export const nodeTests = (dir: string, tier: 'fast' | 'slow' | 'all' = 'fast') => {
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
      tags: [SLOW_TAG],
      tagsFilter: tier === 'fast' ? ['!slow'] : tier === 'slow' ? ['slow'] : undefined,
      projects: mocking.length > 0 ? [shared, isolated] : [shared],
    },
  };
};

/**
 * Heap guards, in Node: a hot path called a million times, heap growth under 1 MB. `heapGrowth` (core/_test/model/heap.ts)
 * fails the run when a GC happened during it; the fixed 64 MB young generation keeps one from hiding allocation.
 */
export const heapTests = () => ({
  ...SOURCES,
  test: {
    include: ['_test/**/*.heap.ts'],
    exclude: [...BROWSER_DIRS, '**/node_modules/**'],
    passWithNoTests: true,
    pool: 'forks',
    fileParallelism: false,
    execArgv: GUARD_FLAGS,
  },
});

/**
 * The quick tier's flags (`test:perf`): a 1 MB young generation (the measurement accounts every collection,
 * core/_test/model/alloc.ts), `--single-threaded` (V8 compiles on the main thread, so tier-up falls on the same step
 * under any load) and lower tier-up thresholds: with V8's defaults a step sits in its middle tier (Maglev) for 15,000
 * to 60,000 steps, allocating there (the wrapped lanes' draw 130 bytes a step, a hover 9 to 26), and a count short enough
 * to be quick ends inside that climb; with these it is over in a few thousand, and the steady state is the defaults'
 * (2026-10-14, Node 24 and 25, M4). `PERF_FLAGS` replaces them.
 */
const PERF_FLAGS = [
  '--expose-gc',
  '--max-semi-space-size=1',
  '--single-threaded',
  '--invocation-count-for-turbofan=300',
  '--minimum-invocations-after-ic-update=50',
];

/**
 * The quick allocation tier (`test:perf`): the `*.quick.ts` cases under `_test` (`quickCase`, core/_test/model/quick.ts),
 * each file in a fork of its own, in parallel and beside any other work: bytes allocated a step, read exactly, against
 * `_test/perf-baseline.json`. Fails a case over its baseline by more than `PERF_THRESHOLD` bytes a step (default 4) or a tenth of the baseline;
 * `--update` (or `PERF_UPDATE=1`) rewrites the baseline. Its own module cache: Vitest clears the shared
 * `node_modules/.vitest-cache` at start, and two runs side by side failed on it (ENOTEMPTY).
 */
export const perfTests = (dir: string) => {
  const baselinePath = join(dir, '_test', 'perf-baseline.json');
  const update = process.env.PERF_UPDATE === '1' || process.argv.includes('--update') || process.argv.includes('-u');
  const threshold = Number(process.env.PERF_THRESHOLD ?? 4);
  return {
    ...SOURCES,
    test: {
      include: ['_test/**/*.quick.ts'],
      exclude: [...BROWSER_DIRS, '**/node_modules/**'],
      passWithNoTests: true,
      pool: 'forks',
      execArgv: process.env.PERF_FLAGS ? process.env.PERF_FLAGS.split(' ') : PERF_FLAGS,
      testTimeout: 120_000,
      hookTimeout: 120_000,
      fsModuleCachePath: join(tmpdir(), `seq-perf-vitest-${process.pid}`),
      reporters: ['default', new PerfReporter(baselinePath, update, threshold)],
      provide: {
        perfBaseline: existsSync(baselinePath) ? JSON.parse(readFileSync(baselinePath, 'utf8')) : {},
        perfUpdate: update,
        perfThreshold: threshold,
      },
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
