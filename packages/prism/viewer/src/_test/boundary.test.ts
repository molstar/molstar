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

// Only `api/` and `sources/` read or change the plugin (operations.md §1, contexts.md §7): elsewhere in the viewer, a
// plugin package is imported for its types only. The mocks (`testing/`), the sequence's first slice (`sequence/`, on
// molstar-poc until plan 2.4) and stories are outside the rule.

const src = path.resolve(import.meta.dirname, '..');
const Allowed = ['api', 'sources', 'testing', 'sequence'];
const Plugin = /^@molstar\/(plugin|plugin-headless|plugin-ui|plugin-extension|mvs)(\/|$)/;

function* sources(dir: string): Generator<string> {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name.startsWith('_')) continue;
      yield* sources(full);
    } else if (/\.tsx?$/.test(entry.name)) yield full;
  }
}

/** Value imports in a file: every `import … from`, `export … from` and `import('…')` that is not type-only. */
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
  for (const m of clean.matchAll(/(?<!typeof\s)import\(\s*['"]([^'"]+)['"]\s*\)/g)) out.push(m[1]);
  return out;
}

describe('the viewer', () => {
  it('reads and changes the plugin only in api/ and sources/', () => {
    const offending: string[] = [];
    for (const file of sources(src)) {
      const rel = path.relative(src, file);
      const top = rel.split(path.sep)[0];
      if (Allowed.includes(top) || rel.includes('.stories.')) continue;
      for (const spec of valueImports(fs.readFileSync(file, 'utf8'))) {
        if (Plugin.test(spec)) offending.push(`${rel} imports ${spec}`);
      }
    }
    expect(offending).toEqual([]);
  });

  it('finds value imports and skips type-only ones', () => {
    expect(valueImports(`import { a } from '@molstar/plugin/context';`)).toEqual(['@molstar/plugin/context']);
    expect(valueImports(`import type { A } from '@molstar/plugin/context';`)).toEqual([]);
    expect(valueImports(`import { type A, type B } from '@molstar/plugin/x';`)).toEqual([]);
    expect(valueImports(`export * from '@molstar/plugin/x';`)).toEqual(['@molstar/plugin/x']);
    expect(valueImports(`const m = await import('@molstar/plugin/x');`)).toEqual(['@molstar/plugin/x']);
  });
});
