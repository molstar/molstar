/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The package's folder graph (packages/seq/AGENTS.md, Imports: "inside a package each folder is an entry"): each folder
// of `src` is an entry (an `index.ts`; `io` and `runtimes` are entries by file, reached as `./io/<file>`), and the
// imports between folders, relative or by the package's own name, never form a cycle. Read from the sources, so a
// planted import fails it.

import { readdirSync, readFileSync } from 'node:fs';
import { dirname, join, relative, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';

const Src = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const Self = '@molstar/seqx-core';
/** Folders whose files are the entries (no index). */
const FileEntries: ReadonlySet<string> = new Set(['io', 'runtimes']);

type Graph = Map<string, Set<string>>;

function sources(dir: string): string[] {
  return readdirSync(dir, { withFileTypes: true }).flatMap((e) =>
    e.name === '_test'
      ? []
      : e.isDirectory()
        ? sources(join(dir, e.name))
        : e.name.endsWith('.ts')
          ? [join(dir, e.name)]
          : [],
  );
}

/** The folder (or top-level file) a specifier names, undefined for an external one. */
function target(spec: string, fromDir: string): string | undefined {
  let rel: string;
  if (spec.startsWith(`${Self}/`)) rel = spec.slice(Self.length + 1);
  else if (spec === Self) return 'index';
  else if (spec.startsWith('.')) rel = relative(Src, resolve(fromDir, spec));
  else return undefined;
  if (rel.startsWith('..')) return undefined;
  return rel.split('/')[0].replace(/\.js$/, '');
}

function folderGraph(files: readonly string[], read: (f: string) => string): Graph {
  const graph: Graph = new Map();
  for (const f of files) {
    const from = relative(Src, f).split('/')[0].replace(/\.ts$/, '');
    for (const m of read(f).matchAll(/(?:from|import)\s*\(?\s*'([^']+)'/g)) {
      const to = target(m[1], dirname(f));
      if (to === undefined || to === '' || to === from) continue;
      let set = graph.get(from);
      if (!set) graph.set(from, (set = new Set()));
      set.add(to);
    }
  }
  return graph;
}

/** One cycle as a path (a -> b -> a), or undefined. */
function findCycle(graph: Graph): string[] | undefined {
  const state = new Map<string, 1 | 2>(); // 1 on the path, 2 done
  const path: string[] = [];
  const visit = (n: string): string[] | undefined => {
    if (state.get(n) === 2) return undefined;
    if (state.get(n) === 1) return [...path.slice(path.indexOf(n)), n];
    state.set(n, 1);
    path.push(n);
    for (const m of graph.get(n) ?? []) {
      const c = visit(m);
      if (c) return c;
    }
    path.pop();
    state.set(n, 2);
    return undefined;
  };
  for (const n of graph.keys()) {
    const c = visit(n);
    if (c) return c;
  }
  return undefined;
}

describe('the core package’s folders', () => {
  const files = sources(Src);
  const folders = readdirSync(Src, { withFileTypes: true })
    .filter((e) => e.isDirectory() && e.name !== '_test')
    .map((e) => e.name);

  it('is each an entry', () => {
    expect(folders.length).toBeGreaterThan(8);
    for (const f of folders) {
      const hasIndex = readdirSync(join(Src, f)).includes('index.ts');
      expect(hasIndex || FileEntries.has(f), `${f} has no index.ts and is not a file-entry folder`).toBe(true);
    }
    for (const f of FileEntries) {
      expect(readdirSync(join(Src, f)).includes('index.ts'), `${f} gained an index: drop it from FileEntries`).toBe(
        false,
      );
    }
  });

  it('import each other without a cycle', () => {
    const graph = folderGraph(files, (f) => readFileSync(f, 'utf8'));
    expect(graph.size).toBeGreaterThan(5);
    expect(findCycle(graph)?.join(' -> ')).toBeUndefined();
  });

  it('finds a planted cycle, by a relative and by a named import', () => {
    const text: Record<string, string> = {
      [join(Src, 'lang/a.ts')]: "import { x } from '../theme/b.js';",
      [join(Src, 'theme/b.ts')]: `import { y } from '${Self}/lang';`,
    };
    const graph = folderGraph(Object.keys(text), (f) => text[f]);
    expect(findCycle(graph)).toEqual(['lang', 'theme', 'lang']);
    // the real graph with one edge planted
    const real = folderGraph(files, (f) => readFileSync(f, 'utf8'));
    real.set('theme', new Set([...(real.get('theme') ?? []), 'lang']));
    expect(findCycle(real)).toBeDefined();
  });
});
