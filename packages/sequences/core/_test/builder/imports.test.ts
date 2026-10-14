/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The builder has no runtime dependencies (B1): from the entry, every module it loads once types are
// erased (Node's own stripping, as the runtime does) is one of its own files; no `@molstar/seq-*` module, nor any other.

import { readFileSync } from 'node:fs';
import { stripTypeScriptTypes } from 'node:module';
import { describe, expect, it } from 'vitest';

const SPECIFIER = /\b(?:import|export)\b[^"';]*?\bfrom\s*["']([^"']+)["']|\bimport\s*\(?\s*["']([^"']+)["']/g;

/** The runtime modules reached from `entry`: own files followed, everything else returned. */
function runtimeImports(entry: URL) {
  const seen = new Set<string>();
  const external = new Set<string>();
  const visit = (url: URL) => {
    if (seen.has(url.href)) return;
    seen.add(url.href);
    const js = stripTypeScriptTypes(readFileSync(url, 'utf8'));
    for (const m of js.matchAll(SPECIFIER)) {
      const spec = m[1] ?? m[2];
      if (spec.startsWith('.')) visit(new URL(spec, url));
      else external.add(spec);
    }
  };
  visit(entry);
  return { files: [...seen].map((f) => f.slice(f.lastIndexOf('/') + 1)), external: [...external] };
}

describe('dependencies', () => {
  it("the entry loads only the builder's own files", () => {
    const { files, external } = runtimeImports(new URL('../../src/builder/index.ts', import.meta.url));
    expect(files.sort()).toEqual(['builder.ts', 'colors.ts', 'helpers.ts', 'index.ts']);
    expect(external).toEqual([]);
  });

  it('the scan sees a value import that stripping keeps, and drops a type-only one', () => {
    const js = (s: string) => [...stripTypeScriptTypes(s).matchAll(SPECIFIER)].map((m) => m[1] ?? m[2]);
    expect(js('import type { A } from "@molstar/seq-core/format";')).toEqual([]);
    expect(js('import { isNamespaced } from "@molstar/seq-core/format";')).toEqual(['@molstar/seq-core/format']);
    expect(js('import { type A, b } from "@molstar/seq-core/model";')).toEqual(['@molstar/seq-core/model']);
  });
});
