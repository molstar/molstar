/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as path from 'node:path';

// Finds the module specifiers a source file loads, for the architecture checks (`dom-free.test.ts`). Text-based, not a
// parser: it knows `import … from`, `export … from`, `export * from`, side-effect `import '…'` and `import('…')`.

export interface ImportRef {
  spec: string;
  /** Erased at compile time: `import type`, `export type … from`, or only `type` names, or `typeof import('…')`. */
  typeOnly: boolean;
}

const SelfPackage = '@molstar/prism-ui';

/** Comments out of the way (a specifier in a comment is not an import); string contents are kept. */
function withoutComments(text: string): string {
  return text.replace(/\/\*[\s\S]*?\*\//g, ' ').replace(/(^|[^:'"`\w])\/\/.*$/gm, '$1');
}

function namesAreTypes(clause: string): boolean {
  const names = clause.match(/\{([^}]*)\}/)?.[1];
  if (names === undefined) return false;
  const list = names
    .split(',')
    .map((n) => n.trim())
    .filter(Boolean);
  return list.length > 0 && list.every((n) => n.startsWith('type '));
}

export function scanImports(source: string): ImportRef[] {
  const text = withoutComments(source);
  const refs: ImportRef[] = [];

  // import … from '…' and export … from '…' (a clause that starts with `*`, `{`, a default name or `type`).
  for (const m of text.matchAll(
    /(?:^|[;{}\s])(import|export)\s+(type\s+)?(\*|\{|[\w$]+)([^;'"`]*?)\s*\bfrom\s*(['"])([^'"\n]+)\5/g,
  )) {
    const [, keyword, typeKeyword, head, rest, , spec] = m;
    // `export const x = a\nfrom` can not match: an export's clause is `*` or `{`.
    if (keyword === 'export' && head !== '*' && head !== '{') continue;
    const clause = head + rest;
    // `import X, { type T }` has a default import: a value. `import { type A }` and `import type …` are erased.
    const typeOnly = !!typeKeyword || (head === '{' && namesAreTypes(clause));
    refs.push({ spec, typeOnly });
  }

  // import '…'
  for (const m of text.matchAll(/(?:^|[;{}\s])import\s*(['"])([^'"\n]+)\1/g))
    refs.push({ spec: m[2], typeOnly: false });

  // import('…'), and `typeof import('…')` in a type position
  for (const m of text.matchAll(/(typeof\s+)?\bimport\s*\(\s*(['"])([^'"\n]+)\2\s*\)/g)) {
    refs.push({ spec: m[3], typeOnly: !!m[1] });
  }
  return refs;
}

/**
 * The top-level folder of `src` a specifier from `file` points into (`dom`, `components`, …; `index` for the package's
 * index): from a relative specifier or a `@molstar/prism-ui/…` one; undefined for any other package.
 */
export function targetFolder(spec: string, file: string, src: string): string | undefined {
  let relative: string;
  if (spec.startsWith('.')) relative = path.relative(src, path.resolve(path.dirname(file), spec));
  else if (spec === SelfPackage) relative = 'index';
  else if (spec.startsWith(`${SelfPackage}/`)) relative = spec.slice(SelfPackage.length + 1);
  else return undefined;
  const top = relative.split(/[\\/]/)[0];
  return top.replace(/\.(?:[cm]?[jt]sx?)$/, '');
}
