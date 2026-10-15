/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import * as fs from 'node:fs';
import * as path from 'node:path';

// What a headless viewer uses is free of the DOM (spec library.md §1, README §3.3): these folders take no value import
// from the folders that touch it (types are fine: a view renderer's type is part of a contract). So no separate core
// package is needed.
const src = path.resolve(import.meta.dirname, '..');
const DomFree = ['reactive', 'commands', 'forms', 'managers', 'model', 'shell'];
const Dom = ['dom', 'react', 'styles', 'components', 'views', 'testing'];

function* sources(dir: string): Generator<string> {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name.startsWith('_')) continue;
      yield* sources(full);
    } else if (/\.tsx?$/.test(entry.name)) yield full;
  }
}

describe('the DOM-free folders', () => {
  it('take no value import from the DOM folders', () => {
    const offending: string[] = [];
    for (const folder of DomFree) {
      for (const file of sources(path.join(src, folder))) {
        if (file.endsWith('.tsx')) offending.push(`${path.relative(src, file)}: a .tsx file`);
        const text = fs.readFileSync(file, 'utf8');
        for (const m of text.matchAll(/^import\s+(type\s+)?[^;]*?from\s+'([^']+)'/gms)) {
          if (m[1]) continue;
          const spec = m[2];
          if (!spec.startsWith('.')) continue;
          const target = path.relative(src, path.resolve(path.dirname(file), spec));
          const top = target.split(path.sep)[0];
          // Type-only named imports (`import { type A, type B }`) are fine too.
          const names = m[0].match(/\{([^}]*)\}/)?.[1];
          const typeOnly = !!names && names.split(',').every((n) => !n.trim() || n.trim().startsWith('type '));
          if (Dom.includes(top) && !typeOnly) offending.push(`${path.relative(src, file)} imports ${spec}`);
        }
      }
    }
    expect(offending).toEqual([]);
  });
});
