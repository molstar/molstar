/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { globSync, readFileSync } from 'node:fs';
import { basename } from 'node:path';
import { defaultClientConditions, defaultServerConditions } from 'vite';
import { defineConfig } from 'vitest/config';

const include = ['**/*.test.{ts,tsx,mjs}'];
// as `exclude`: these folders at any depth, these at the root
const skipDirs = new Set(['node_modules', 'lib', 'build']);
const skipRoots = new Set(['.pnpm-store', '.claude', '.git', 'smoke']);
const exclude = ['**/node_modules/**', '**/lib/**', '**/build/**', '.pnpm-store/**', '.claude/**', 'smoke/**'];

// The test lanes (AGENTS.md, Tests), by Vitest tag: a test is quick unless it carries `{ tags: ['medium'] }` or
// `{ tags: ['heavy'] }` (a `describe` tags its tests), or its file `// @module-tag medium`. `vitest run` is the quick
// lane, `--mode medium` and `--mode heavy` the others, `--mode pr` quick and medium (what a PR's checks run), `--mode
// all` every lane. Tags inherit: a heavy test inside a medium `describe` is heavy.
const Lanes = {
  quick: ['!medium', '!heavy'],
  medium: ['medium', '!heavy'],
  heavy: ['heavy'],
  pr: ['!heavy'],
  all: undefined,
} satisfies Record<string, string[] | undefined>;
type Lane = keyof typeof Lanes;

const tags = [
  { name: 'medium', description: 'larger fixtures, the native GL suite, sweeps: `pnpm test:medium`, its own CI job' },
  {
    name: 'heavy',
    description: 'wall-clock and scaling bounds, full-size fuzz, CLIs, Chrome: `pnpm test:heavy`, nightly or on demand',
    timeout: 120_000,
  },
];

// Every test file, read once: its module tag, whether it names the medium tag, whether it mocks modules (Seq*). A
// lane only imports the files that can hold its tests (a filtered-out file is still imported, and imports are about
// half the run): the quick lane skips the files tagged medium or heavy as a whole, the medium lane takes the files that
// name the tag in an array (`tags: ['medium']`: medium tags sit in the test's own file, never in a helper; a test file
// imported by another runs its medium tests in its own file only, as `points.test.ts` in `renderer.test.ts`). The heavy
// lane imports every file: helpers tag its tests (`itWithin`, Seq*'s smoke stories), and it runs nightly.
const Mocks = /\bvi\.(mock|doMock|unmock|resetModules)\(/;
const files = globSync(include, {
  exclude: (p: string) => skipDirs.has(basename(p)) || skipRoots.has(p.replaceAll('\\', '/')),
})
  .map((f) => f.replaceAll('\\', '/'))
  .sort()
  .map((path) => {
    const source = readFileSync(path, 'utf8');
    return {
      path,
      module: /^\s*(?:\/\/|\/?\*+)\s*@module-tag\s+(medium|heavy)\b/m.exec(source)?.[1],
      medium: /\[\s*['"]medium['"]/.test(source),
      mocks: Mocks.test(source),
    };
  });

// Mol*'s files run isolated, each in a module registry of its own. Shared registries (`isolate: false`, as Seq* runs)
// would save about a third of the quick lane's CPU (Mol*'s quick files, M4, 2026-10-20: 672 CPU-s isolated, 423
// shared), but which file a worker ran first then decides how an import cycle evaluates: under Vite's module runner an
// `export *` copies the names the module has at that moment, so `to_mmCIF` (`model/structure`, re-exported from
// `export/mmcif.js`) is undefined after `formats/structure/mmcif` was imported first; and module state (registries,
// spies on namespaces) leaks between files. Sharding the CI job buys the same time without the flakes.

// Seq* (`packages/sequences`, its AGENTS.md, Tests) runs here in two projects of its own, under its own rules: files
// share a worker's module registry (`isolate: false`), the files that mock or reset modules run isolated. Never
// matched: the browser folders, the heap guards (`*.heap.ts`, not `*.test.*`) and the benches; its opt-in tiers
// (`test:perf`, `test:browser`, `test:heap`, `bench`) are the packages' own scripts (`sequences/vitest.shared.ts`).
const seq = 'packages/sequences/';
const seqBrowser = /^packages\/sequences\/[^/]+\/(_test\/(.*\/)?(canvas|dom|react)\/|bench\/)/;

function lane(mode: string): Lane {
  if (mode === 'test') return 'quick';
  if (mode in Lanes) return mode as Lane;
  throw new Error(`Unknown test lane '${mode}': use one of ${Object.keys(Lanes).join(', ')} (vitest --mode <lane>).`);
}

function inLane(l: Lane, f: (typeof files)[number]) {
  if (l === 'quick') return !f.module;
  if (l === 'medium') return f.module === 'medium' || f.medium;
  if (l === 'pr') return f.module !== 'heavy';
  return true;
}

// One Vitest run for the whole workspace, ESM throughout. Workspace packages resolve to their sources through the
// `molstar-src` export condition (the same one the build and the type checker use), so a new package needs no mapping
// here. DOM tests opt in per file with `// @vitest-environment jsdom`.
export default defineConfig(({ mode }) => {
  const l = lane(mode);
  const selected = files.filter((f) => inLane(l, f) && !seqBrowser.test(f.path));
  const project = (name: string, list: string[], isolate: boolean) => ({
    extends: true as const,
    test: { name, include: list, exclude, isolate, sequence: { groupOrder: 0 } },
  });
  const seqFiles = selected.filter((f) => f.path.startsWith(seq));
  return {
    resolve: { conditions: ['molstar-src', ...defaultClientConditions] },
    ssr: { resolve: { conditions: ['molstar-src', ...defaultServerConditions] } },
    test: {
      // Preact (the adapter tests) ships ESM and CommonJS builds; inlined, all of it comes from one.
      server: { deps: { inline: ['preact'] } },
      // The long cases (2,000-append histories, fresh module graphs) can share the machine with other files.
      testTimeout: 30_000,
      hookTimeout: 30_000,
      environmentOptions: { jsdom: { url: 'http://localhost/' } },
      tags,
      tagsFilter: Lanes[l],
      // a filtered run of the medium or heavy lane may hold none of their tests; the quick lane always has some
      passWithNoTests: l === 'medium' || l === 'heavy',
      projects: [
        project(
          'unit',
          selected.filter((f) => !f.path.startsWith(seq)).map((f) => f.path),
          true,
        ),
        // (a file that mocks or resets modules needs a module registry of its own)
        ...[
          project(
            'sequences',
            seqFiles.filter((f) => !f.mocks).map((f) => f.path),
            false,
          ),
          project(
            'sequences-isolated',
            seqFiles.filter((f) => f.mocks).map((f) => f.path),
            true,
          ),
        ].filter((p) => p.test.include.length > 0),
      ],
    },
  };
});
