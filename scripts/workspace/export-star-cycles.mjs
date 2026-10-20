/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// No `export * from` inside an import cycle. Native ESM binds a star export's names at link time, but Vite's module
// runner (Vitest, the dev server) copies the names the target has when the `export *` runs: inside a cycle that can
// be before the target ran, and the names are then missing for good, depending on which module of the cycle was
// imported first (`to_mmCIF` from `model/structure` was undefined after `formats/structure/mmcif`). Tests that share a
// module registry (AGENTS.md, Tests) would pass or fail by the order of their files. Name the exports instead:
// `export { a, type B } from './x.js'` reads them live. Fast enough for `pnpm lint`.
//
// Usage: node scripts/workspace/export-star-cycles.mjs

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { buildImportGraph } from './import-graph.mjs';

/** The strongly connected components of a graph (Tarjan, iterative), as a map from each node to its component. */
export function components(graph) {
  const index = new Map();
  const low = new Map();
  const onStack = new Set();
  const stack = [];
  const component = new Map();
  let next = 0;
  let count = 0;
  const edges = (v) => graph.get(v) ?? new Set();
  for (const root of graph.keys()) {
    if (index.has(root)) continue;
    const work = [[root, edges(root)[Symbol.iterator]()]];
    index.set(root, next);
    low.set(root, next++);
    stack.push(root);
    onStack.add(root);
    while (work.length) {
      const [v, it] = work[work.length - 1];
      const step = it.next();
      if (!step.done) {
        const w = step.value;
        if (!index.has(w)) {
          index.set(w, next);
          low.set(w, next++);
          stack.push(w);
          onStack.add(w);
          work.push([w, edges(w)[Symbol.iterator]()]);
        } else if (onStack.has(w)) low.set(v, Math.min(low.get(v), index.get(w)));
        continue;
      }
      work.pop();
      if (work.length) {
        const parent = work[work.length - 1][0];
        low.set(parent, Math.min(low.get(parent), low.get(v)));
      }
      if (low.get(v) === index.get(v)) {
        let w;
        do {
          w = stack.pop();
          onStack.delete(w);
          component.set(w, count);
        } while (w !== v);
        count++;
      }
    }
  }
  return component;
}

/** The `export *` edges whose two ends lie on one import cycle (evaluated imports only: no `import()`). */
export function starsInCycles(graph) {
  const component = components(graph.static);
  return graph.stars.filter(({ from, to }) => component.has(from) && component.get(from) === component.get(to));
}

const isMain = process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url);
if (isMain) {
  const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
  const inventory = JSON.parse(fs.readFileSync(path.join(root, 'scripts/workspace/inventory.json'), 'utf8'));
  const found = starsInCycles(buildImportGraph({ root, packages: inventory.packages ?? [] }));
  if (found.length) {
    console.error('`export * from` inside an import cycle (scripts/workspace/export-star-cycles.mjs says why):');
    for (const { from, to } of found) console.error(`  ${from} -> ${to}`);
    console.error("Name the exports instead: export { a, type B } from './x.js'.");
    process.exitCode = 1;
  }
}
