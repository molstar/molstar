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
import { scanImports, targetFolder } from './import-scan.js';

// What a headless viewer uses is free of the DOM (spec library.md §1, api.md §11): these folders take no value import
// from the folders that touch it (types are fine: a view renderer's type is part of a contract). So no separate core
// package is needed.
const src = path.resolve(import.meta.dirname, '..');
const DomFree = ['reactive', 'commands', 'forms', 'managers', 'model', 'shell'];
// `index` is the package's index, which re-exports the DOM folders.
const Dom = ['dom', 'react', 'styles', 'components', 'views', 'testing', 'index'];

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
        const name = path.relative(src, file);
        if (file.endsWith('.tsx')) offending.push(`${name}: a .tsx file`);
        for (const ref of scanImports(fs.readFileSync(file, 'utf8'))) {
          const top = targetFolder(ref.spec, file, src);
          if (top && Dom.includes(top) && !ref.typeOnly) offending.push(`${name} imports ${ref.spec}`);
        }
      }
    }
    expect(offending).toEqual([]);
  });
});

describe('the import scan', () => {
  const scan = (text: string) => scanImports(text).map((r) => `${r.typeOnly ? 'type ' : ''}${r.spec}`);

  it('finds import and export forms', () => {
    expect(scan(`import { a } from '../dom/x.js';`)).toEqual(['../dom/x.js']);
    expect(scan(`import X from '../dom/x.js'`)).toEqual(['../dom/x.js']);
    expect(scan(`import * as ns from "../dom/x.js";`)).toEqual(['../dom/x.js']);
    expect(scan(`export { a, b } from '../dom/x.js';`)).toEqual(['../dom/x.js']);
    expect(scan(`export * from '../dom/x.js';`)).toEqual(['../dom/x.js']);
    expect(scan(`export * as ns from '../dom/x.js';`)).toEqual(['../dom/x.js']);
    expect(scan(`import '../dom/register.js';`)).toEqual(['../dom/register.js']);
    expect(scan(`const m = await import('../dom/x.js');`)).toEqual(['../dom/x.js']);
  });

  it('spans lines', () => {
    expect(
      scan(`import {\n  a,\n  b,\n} from '../dom/x.js';\nimport {\n  type A,\n  type B,\n} from '../dom/y.js';`),
    ).toEqual(['../dom/x.js', 'type ../dom/y.js']);
  });

  it('knows what the compiler erases', () => {
    expect(scan(`import type { A } from '../dom/x.js';`)).toEqual(['type ../dom/x.js']);
    expect(scan(`import type X from '../dom/x.js';`)).toEqual(['type ../dom/x.js']);
    expect(scan(`import type * as ns from '../dom/x.js';`)).toEqual(['type ../dom/x.js']);
    expect(scan(`import { type A, type B } from '../dom/x.js';`)).toEqual(['type ../dom/x.js']);
    expect(scan(`export type { A } from '../dom/x.js';`)).toEqual(['type ../dom/x.js']);
    expect(scan(`export { type A } from '../dom/x.js';`)).toEqual(['type ../dom/x.js']);
    expect(scan(`type T = typeof import('../dom/x.js');`)).toEqual(['type ../dom/x.js']);
  });

  it('keeps mixed imports as values', () => {
    expect(scan(`import X, { type T } from '../dom/x.js';`)).toEqual(['../dom/x.js']);
    expect(scan(`import { type A, b } from '../dom/x.js';`)).toEqual(['../dom/x.js']);
    expect(scan(`export { type A, b } from '../dom/x.js';`)).toEqual(['../dom/x.js']);
  });

  it('ignores comments and code that only looks like an import', () => {
    expect(scan(`// import { a } from '../dom/x.js';\n/**\n * import '../dom/y.js'\n */`)).toEqual([]);
    expect(scan(`export const from = 1;\nconst s = "import x from 'y'";`)).toEqual([]);
    expect(scan(`const url = 'http://a/b'; import { a } from '../x.js';`)).toEqual(['../x.js']);
  });

  it('resolves the target folder from relative and self-package specifiers', () => {
    const file = path.join(src, 'managers', 'a', 'b.ts');
    expect(targetFolder('../../dom/x.js', file, src)).toBe('dom');
    expect(targetFolder('./c.js', file, src)).toBe('managers');
    expect(targetFolder('../../index.js', file, src)).toBe('index');
    expect(targetFolder('@molstar/prism-ui/components/button/button.js', file, src)).toBe('components');
    expect(targetFolder('@molstar/prism-ui/dom', file, src)).toBe('dom');
    expect(targetFolder('@molstar/prism-ui/dom.js', file, src)).toBe('dom');
    expect(targetFolder('@molstar/prism-ui', file, src)).toBe('index');
    expect(targetFolder('rxjs', file, src)).toBeUndefined();
    expect(targetFolder('@molstar/prism-ui-other/dom', file, src)).toBeUndefined();
  });
});
