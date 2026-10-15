/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's edge: nothing in src/ reaches outside it by a path, and its only dependencies are the core's entries,
// Mol*'s core (its input and binding vocabulary), rxjs, React (the /react entry), Node's built-ins and Playwright (the
// CLI's rasterizer, loaded when present). Its inner layering is layers.grit's.

import { readdirSync, readFileSync, statSync } from 'node:fs';
import { dirname, join, relative, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';

const SRC = fileURLToPath(new URL('../src', import.meta.url));
const EXTERNAL =
  /^(?:@molstar\/seq-core\/|@molstar\/seq-viewer(?:\/|$)|@molstar\/core\/|rxjs(?:\/|$)|react(?:-dom)?(?:\/|$)|node:|playwright$)/;
/** Static imports and re-exports (up to their `from`), bare imports, and dynamic imports. */
const SPECIFIER =
  /^(?:import|export)\s[^;]*?[\s}*]from\s*["']([^"']+)["']|^import\s*["']([^"']+)["']|\bimport\(\s*["']([^"']+)["']/gm;

function files(dir: string, out: string[] = []): string[] {
  for (const name of readdirSync(dir)) {
    const f = join(dir, name);
    if (statSync(f).isDirectory()) files(f, out);
    else if (/\.tsx?$/.test(name)) out.push(f);
  }
  return out;
}

describe("the viewer's edge", () => {
  it('src/ imports by path only inside itself, and only the dependencies it declares', () => {
    const all = files(SRC);
    expect(all.length).toBeGreaterThan(150);
    const wrong: string[] = [];
    for (const file of all) {
      for (const m of readFileSync(file, 'utf8').matchAll(SPECIFIER)) {
        const spec = m[1] ?? m[2] ?? m[3];
        const at = relative(SRC, file);
        if (spec.startsWith('.')) {
          if (relative(SRC, resolve(dirname(file), spec)).startsWith('..')) wrong.push(`${at}: ${spec} leaves src/`);
        } else if (!EXTERNAL.test(spec)) wrong.push(`${at}: ${spec}`);
      }
    }
    expect(wrong).toEqual([]);
  });
});
