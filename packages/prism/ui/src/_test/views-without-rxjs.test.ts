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

// Views do no rxjs (rendering.md §2.6): they get states and derived values and render them (`watch`, `bind` by
// `listen`); rxjs is the models' and sources' boundary (`toState`, `asObservable`). Strict (2026-10-17, the views all
// migrated: 357 uses in 42 files at the start): an rxjs import (type-only too) or a `.pipe(` in a view fails.

const prism = path.resolve(import.meta.dirname, '../../..');
const repo = path.resolve(prism, '../..');

function* files(dir: string): Generator<string> {
  if (!fs.existsSync(dir)) return;
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name === '_test' || entry.name === '_stories' || entry.name === '_generated') continue;
      yield* files(full);
    } else if (/\.tsx?$/.test(entry.name) && !entry.name.endsWith('.d.ts')) yield full;
  }
}

/** The view files: views/ in ui, viewer and assistant, ui's components, the extensions' Prism UI (.tsx on prism-ui). */
function viewFiles(): string[] {
  const out: string[] = [];
  const viewsUnder = (dir: string) => {
    for (const f of files(dir)) if (f.split(path.sep).includes('views')) out.push(f);
  };
  for (const pkg of ['ui', 'viewer', 'assistant']) viewsUnder(path.join(prism, pkg, 'src'));
  for (const f of files(path.join(prism, 'ui/src/components'))) out.push(f);
  const extensions = path.join(repo, 'extensions');
  for (const ext of fs.readdirSync(extensions)) {
    for (const f of files(path.join(extensions, ext, 'src'))) {
      if (f.endsWith('.tsx') && fs.readFileSync(f, 'utf8').includes('@molstar/prism-ui')) out.push(f);
    }
  }
  return out;
}

/** rxjs imports (type-only ones too) and `.pipe(` calls in a file. */
export function rxjsUse(source: string): number {
  const text = source.replace(/(^|[^:])\/\/.*$/gm, '$1');
  const imports = [...text.matchAll(/from\s+['"]rxjs(\/[^'"]*)?['"]/g)].length;
  const pipes = [...text.matchAll(/\.pipe\(/g)].length;
  return imports + pipes;
}

describe('views without rxjs', () => {
  it('counts imports and pipes', () => {
    expect(rxjsUse("import { map } from 'rxjs';\nx.pipe(map(f));\n// y.pipe(\n")).toBe(2);
    expect(rxjsUse("import { derive } from '@molstar/prism-ui/reactive/state';\nderive(x, f);")).toBe(0);
  });

  it('no rxjs in a view', () => {
    const found: string[] = [];
    for (const f of viewFiles()) {
      const n = rxjsUse(fs.readFileSync(f, 'utf8'));
      if (n > 0) found.push(`${path.relative(repo, f).split(path.sep).join('/')}: ${n}`);
    }
    expect(found, 'rxjs in a view: use states, derive and watch (rendering.md §2.6)').toEqual([]);
  });
});
