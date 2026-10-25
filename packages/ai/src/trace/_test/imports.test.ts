/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readFileSync } from 'node:fs';
import { dirname, join, relative, sep } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';

const src = join(dirname(fileURLToPath(import.meta.url)), '../..');

/** Every source file a module reaches through relative imports, as paths under `src/`. */
function reached(entry: string): Set<string> {
  const seen = new Set<string>();
  const visit = (file: string) => {
    const key = relative(src, file).split(sep).join('/');
    if (seen.has(key)) return;
    seen.add(key);
    const text = readFileSync(file, 'utf8');
    for (const [, spec] of text.matchAll(/(?:import|export)[^'"]*?from\s+'(\.[^']+)'/g))
      visit(join(dirname(file), (spec as string).replace(/\.js$/, '.ts')));
  };
  visit(join(src, entry));
  return seen;
}

describe('the trace entry', () => {
  it('reaches neither the session nor the tools (so not zod): the markers are in messages/', () => {
    const files = [...reached('trace/index.ts')];
    expect(files.filter((f) => f.startsWith('session/') || f.startsWith('tools/'))).toEqual([]);
    expect(files).toContain('messages/markers.ts');
  });
});
