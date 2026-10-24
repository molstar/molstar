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

// Every test file, read once: its module tag, whether it names the medium tag, whether it mocks modules. A
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

// Files share their worker's module registry (`isolate: false`): each file importing the libraries again was most of
// the run (Mol*'s quick files, M4, 2026-10-20: 672 CPU-s isolated, 423 shared; Seq* found the same,
// packages/sequences/AGENTS.md). Two rules keep a file's result independent of the files before it in its worker:
// no `export *` inside an import cycle (`scripts/workspace/export-star-cycles.mjs`, in `pnpm lint`), and module state a
// file changes it restores (`vitest.setup.ts` restores spies, stubbed globals and environment variables and real timers
// after every file). A file that mocks or resets modules runs isolated, found by reading it; so do the files below,
// which need a registry of their own (shuffled runs find more: `vitest run --sequence.shuffle.files --maxWorkers=3`).
const Isolated = new Set([
  // fake frame timers must not inherit a pending tick from another file's discarded clock
  'packages/prism/ui/src/views/canvas/_test/hover-label-frames.test.tsx',
  // asserts which languages are registered: another file's import registers them too
  'packages/query-language/src/_test/compile.test.ts',
  // lists every registered provider: another file's registrations are listed too
  'cli/state-docs/src/_test/state-docs.test.ts',
  // spies on `Vec3.orthogonalize` before `curve-segment.ts` loads and keeps it in a constant
  'packages/graphics/src/repr/structure/visual/util/polymer/_test/cartoon-geometry.test.ts',
  // registers a stand-in encoder in Mediabunny's registry, which has no unregister
  'packages/ai/src/narration/_test/mediabunny-containers.test.ts',
]);
const isolated = (f: (typeof files)[number]) => f.mocks || Isolated.has(f.path);
const setupFiles = [`${import.meta.dirname}/vitest.setup.ts`];

// Seq* (`packages/sequences`, its AGENTS.md, Tests) runs here in two projects of its own, under the same rules. Never
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
  const project = (name: string, list: typeof selected, isolate: boolean) => ({
    extends: true as const,
    test: {
      name,
      include: list.map((f) => f.path),
      exclude,
      isolate,
      ...(isolate ? {} : { setupFiles }),
      sequence: { groupOrder: 0 },
    },
  });
  const molstar = selected.filter((f) => !f.path.startsWith(seq));
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
          molstar.filter((f) => !isolated(f)),
          false,
        ),
        ...[
          project('unit-isolated', molstar.filter(isolated), true),
          project(
            'sequences',
            seqFiles.filter((f) => !isolated(f)),
            false,
          ),
          project('sequences-isolated', seqFiles.filter(isolated), true),
        ].filter((p) => p.test.include.length > 0),
      ],
    },
  };
});
