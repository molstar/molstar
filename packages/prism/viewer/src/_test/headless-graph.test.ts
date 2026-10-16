/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { describe, expect, it } from 'vitest';

// The headless viewer's module graph is free of the DOM (api.md §11, library.md §1): following its value imports from
// `headless.ts` through the viewer and into prism-ui never reaches prism-ui's DOM folders or its index (which
// re-exports them); prism-ui's DOM-free part is `@molstar/prism-ui/core`.

const viewerSrc = path.resolve(import.meta.dirname, '..');
const uiSrc = path.resolve(import.meta.dirname, '../../../ui/src');
const UiDom = ['dom', 'react', 'styles', 'components', 'views', 'testing'];

/** Value import specifiers (type-only ones are erased by the compiler). */
function valueImports(text: string): string[] {
  const out: string[] = [];
  const clean = text.replace(/\/\*[\s\S]*?\*\//g, '').replace(/^\s*\/\/.*$/gm, '');
  for (const m of clean.matchAll(/(?:^|;|\n)\s*(import|export)\s+(type\s+)?([\s\S]*?)\s+from\s+['"]([^'"]+)['"]/g)) {
    if (m[2]) continue;
    const names = m[3].match(/^\{([^}]*)\}$/)?.[1];
    if (names && names.split(',').every((n) => !n.trim() || n.trim().startsWith('type '))) continue;
    out.push(m[4]);
  }
  for (const m of clean.matchAll(/(?:^|;|\n)\s*import\s+['"]([^'"]+)['"]/g)) out.push(m[1]);
  return out;
}

function resolveFile(spec: string, from: string): string | undefined {
  let base: string | undefined;
  if (spec.startsWith('.')) base = path.resolve(path.dirname(from), spec);
  else if (spec === '@molstar/prism-ui') base = path.join(uiSrc, 'index.js');
  else if (spec.startsWith('@molstar/prism-ui/')) base = path.join(uiSrc, spec.slice('@molstar/prism-ui/'.length));
  else if (spec.startsWith('@molstar/prism-viewer/'))
    base = path.join(viewerSrc, spec.slice('@molstar/prism-viewer/'.length));
  if (!base) return undefined;
  const stem = base.replace(/\.js$/, '');
  for (const candidate of [`${stem}.ts`, `${stem}.tsx`, path.join(stem, 'index.ts')]) {
    if (fs.existsSync(candidate)) return candidate;
  }
  return undefined;
}

describe('the headless viewer', () => {
  it('reaches none of prism-ui DOM modules by value', () => {
    const seen = new Set<string>();
    const offending: string[] = [];
    const queue: [string, string[]][] = [[path.join(viewerSrc, 'headless.ts'), []]];
    while (queue.length) {
      const [file, chain] = queue.shift()!;
      if (seen.has(file)) continue;
      seen.add(file);
      const relUi = path.relative(uiSrc, file);
      if (!relUi.startsWith('..')) {
        const top = relUi.split(path.sep)[0];
        if (UiDom.includes(top) || relUi === 'index.ts') {
          offending.push([...chain, path.relative(path.resolve(uiSrc, '..', '..'), file)].join(' -> '));
          continue;
        }
      }
      for (const spec of valueImports(fs.readFileSync(file, 'utf8'))) {
        const next = resolveFile(spec, file);
        if (next) queue.push([next, [...chain, path.relative(path.resolve(uiSrc, '..', '..'), file)]]);
      }
    }
    expect(offending.slice(0, 5)).toEqual([]);
    expect(seen.size).toBeGreaterThan(20);
  });
});
