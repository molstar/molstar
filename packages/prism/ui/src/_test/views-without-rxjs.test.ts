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
import { RxjsInViews } from './views-without-rxjs.offenders.js';

// Views do no rxjs (rendering.md §2.6): they get states and derived values and render them (`watch`, `bind` by
// `listen`); rxjs is the models' and sources' boundary (`toState`, `asObservable`). Phased in: the files that still
// import rxjs or call `.pipe(` are listed with their counts (views-without-rxjs.offenders.ts); a new one, or more in a
// listed one, fails; a listed one that has fewer must lower its entry (the list only shrinks).

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

  it('no new rxjs in views; the listed ones only shrink', () => {
    const found: Record<string, number> = {};
    for (const f of viewFiles()) {
      const n = rxjsUse(fs.readFileSync(f, 'utf8'));
      if (n > 0) found[path.relative(repo, f).split(path.sep).join('/')] = n;
    }
    const grown: string[] = [];
    const shrunk: string[] = [];
    for (const [file, n] of Object.entries(found)) {
      const allowed = RxjsInViews[file] ?? 0;
      if (n > allowed) grown.push(`${file}: ${n} (allowed ${allowed})`);
      else if (n < allowed) shrunk.push(`${file}: ${n} (listed ${allowed})`);
    }
    for (const [file, n] of Object.entries(RxjsInViews)) if (!(file in found)) shrunk.push(`${file}: 0 (listed ${n})`);
    expect(grown, 'rxjs in a view: use states, derive and watch (rendering.md §2.6)').toEqual([]);
    expect(shrunk, 'fewer than listed: lower or remove the entry in views-without-rxjs.offenders.ts').toEqual([]);
  });
});
