/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The vitest configs every package shares (RB0). A package's own configs call these, so `mise run test` runs the same
// kinds in each: the Node suite (the fast tier, `test`, the one CI runs), the slow Node tier (`test:slow`), the heap
// guards, the browser suite (and its heap guards, `HEAP=1`) and the bench. They return plain objects: this folder has
// no vitest of its own to import from.

import { globSync, readFileSync } from "node:fs";
import { join } from "node:path";

const GUARD_FLAGS = ["--expose-gc", "--min-semi-space-size=64", "--max-semi-space-size=64"];
const BROWSER_DIRS = ["_test/**/canvas/**", "_test/**/dom/**", "_test/**/react/**"];
/** A file that mocks or resets modules needs a module registry of its own (`isolate`). */
const MOCKS = /\bvi\.(mock|doMock|unmock|resetModules)\(/;

/**
 * The slow tier's tag. A test is slow by nature (a wall-clock timing, a fuzz or oracle run of thousands, a million
 * rows, a large fixture, Chrome or a CLI in a process of its own) when it carries `{ tags: ["slow"] }`, or its file
 * `// @module-tag slow`; its quick representative (a short fuzz, a few stories, a small scale) stays untagged.
 */
const SLOW_TAG = {
  name: "slow",
  description: "slow by nature (timing, fuzz, scale, Chrome, CLIs): `test:slow`, never CI",
  timeout: 120_000,
};

/**
 * The Node suites, with no DOM (browser tests run in `browserTests`). `tier`: `"fast"` is `test`, the default and the
 * one CI runs: fast, deterministic, no wall-clock timing, no Chrome, no heap guards; `"slow"` is `test:slow`, the tests
 * tagged slow, run by hand or at checkpoints; `"all"` both. Files share their worker's module registry (`isolate:
 * false`: each file imported the whole package again, most of the run's time); the files that mock or reset modules,
 * found by reading them, run isolated in a project of their own. `dir` is the package's folder.
 */
export const nodeTests = (dir: string, tier: "fast" | "slow" | "all" = "fast") => {
  const exclude = [...BROWSER_DIRS, "**/node_modules/**"];
  const mocking = globSync("_test/**/*.test.ts", { cwd: dir, exclude })
    .filter((f) => MOCKS.test(readFileSync(join(dir, f), "utf8")))
    .sort();
  // A project extending the root concatenates its lists with the root's: the root names no files of its own. Vitest
  // queues projects by name: "isolated" first, so its fresh workers start with the run, not after it.
  const shared = {
    extends: true,
    test: { name: "shared", include: ["_test/**/*.test.ts"], exclude: [...exclude, ...mocking], isolate: false },
  };
  const isolated = { extends: true, test: { name: "isolated", include: mocking } };
  return {
    test: {
      passWithNoTests: true,
      ...(process.env.PROF_SETUP ? { setupFiles: [process.env.PROF_SETUP] } : {}),
      tags: [SLOW_TAG],
      tagsFilter: tier === "fast" ? ["!slow"] : tier === "slow" ? ["slow"] : undefined,
      projects: mocking.length > 0 ? [shared, isolated] : [shared],
    },
  };
};

/**
 * Heap guards, in Node: a hot path called a million times, heap growth under 1 MB. `heapGrowth` (core/_test/model/heap.ts)
 * fails the run when a GC happened during it; the fixed 64 MB young generation keeps one from hiding allocation.
 */
export const heapTests = () => ({
  test: {
    include: ["_test/**/*.heap.ts"],
    exclude: [...BROWSER_DIRS, "**/node_modules/**"],
    passWithNoTests: true,
    pool: "forks",
    fileParallelism: false,
    execArgv: GUARD_FLAGS,
  },
});

/**
 * Chromium on the installed Chrome (no browser download). `HEAP=1` runs the `*.heap.ts` guards instead, with `gc`
 * exposed, the fixed young generation and precise `performance.memory` (_test/canvas/heap-browser.ts).
 * `playwright` is the package's own `@vitest/browser-playwright` export.
 */
export const browserTests = (playwright: (o: object) => unknown) => {
  const heap = process.env.HEAP === "1";
  return {
    test: {
      include: heap ? ["_test/**/{canvas,dom,react}/**/*.heap.ts"] : ["_test/**/{canvas,dom,react}/**/*.test.{ts,tsx}"],
      passWithNoTests: true,
      fileParallelism: false,
      browser: {
        enabled: true,
        headless: true,
        screenshotFailures: false,
        provider: playwright({
          launchOptions: {
            channel: "chrome",
            args: heap ? [`--js-flags=${GUARD_FLAGS.join(" ")}`, "--enable-precise-memory-info"] : [],
          },
        }),
        instances: [{ browser: "chromium", viewport: { width: 1200, height: 900 } }],
      },
    },
  };
};

/** `vitest bench`, outside `test`: the perf phases' numbers (`pnpm --filter <package> bench`). */
export const benchTests = () => ({ test: { benchmark: { include: ["bench/**/*.bench.ts"] } } });
