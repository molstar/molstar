/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The vitest configs of the opt-in guard tiers every package shares (the root AGENTS.md, Tests): the heap guards
// (`heapTests`) and the quick allocation tier (`perfTests`). A package's own `vitest.heap.config.ts` and
// `vitest.perf.config.ts` call these. They return plain objects, not `defineConfig`'s.

import { existsSync, readFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { defaultClientConditions, defaultServerConditions } from 'vite';
import { PerfReporter } from './perf-reporter.js';

/** The heap guards' flags: `gc` exposed, a fixed 64 MB young generation (a collection would hide allocation) */
export const GUARD_FLAGS = ['--expose-gc', '--min-semi-space-size=64', '--max-semi-space-size=64'];
/** The test folders a browser runs (Seq*'s `browserTests`), left out of the Node tiers */
export const BROWSER_DIRS = ['_test/**/canvas/**', '_test/**/dom/**', '_test/**/react/**'];

/** Workspace packages resolve to their sources through the `molstar-src` export condition, as in the root config. */
export const SOURCES = {
  resolve: { conditions: ['molstar-src', ...defaultClientConditions] },
  ssr: { resolve: { conditions: ['molstar-src', ...defaultServerConditions] } },
};

/**
 * Heap guards, in Node: a hot path called a million times, heap growth under 1 MB. `heapGrowth` (heap.ts) fails the run
 * when a GC happened during it; the fixed 64 MB young generation keeps one from hiding allocation.
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
 * The quick tier's flags (`test:perf`): a 1 MB young generation (the measurement accounts every collection, alloc.ts),
 * `--single-threaded` (V8 compiles on the main thread, so tier-up falls on the same step under any load) and lower
 * tier-up thresholds: with V8's defaults a step sits in its middle tier (Maglev) for 15,000 to 60,000 steps, allocating
 * there (the wrapped lanes' draw 130 bytes a step, a hover 9 to 26), and a count short enough to be quick ends inside
 * that climb; with these it is over in a few thousand, and the steady state is the defaults' (2026-10-14, Node 24 and
 * 25, M4). `PERF_FLAGS` replaces them.
 */
const PERF_FLAGS = [
  '--expose-gc',
  '--max-semi-space-size=1',
  '--single-threaded',
  '--invocation-count-for-turbofan=300',
  '--minimum-invocations-after-ic-update=50',
];

/**
 * The quick allocation tier (`test:perf`): the `*.quick.ts` cases under `_test` (`quickCase`, quick.ts), each file in a
 * fork of its own, in parallel and beside any other work: bytes allocated a step, read exactly, against
 * `_test/perf-baseline.json`. Fails a case over its baseline by more than `PERF_THRESHOLD` bytes a step (default 4) or a
 * tenth of the baseline; `--update` (or `PERF_UPDATE=1`) rewrites the baseline. Its own module cache: Vitest clears the
 * shared `node_modules/.vitest-cache` at start, and two runs side by side failed on it (ENOTEMPTY).
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
