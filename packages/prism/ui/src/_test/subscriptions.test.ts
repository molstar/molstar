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

// A component follows a stream with `bind` (dom/scope.ts), which its render scope unsubscribes. A raw `.subscribe(`
// in a component outlives it unless it is kept in a const and unsubscribed (`onCleanup(() => sub.unsubscribe())`):
// the command menu's rows and selected rows were, and a closed menu kept recomputing on the registry's changes.
const src = path.resolve(import.meta.dirname, '..');

function* components(dir: string): Generator<string> {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name === '_test' || entry.name === '_stories' || entry.name === 'testing') continue;
      yield* components(full);
    } else if (entry.name.endsWith('.tsx')) yield full;
  }
}

/** The raw subscriptions in a file that are not kept in a const it unsubscribes. */
function unowned(source: string): string[] {
  // Line comments out (not `://`), so their words and semicolons are not read as code.
  const text = source.replace(/(^|[^:])\/\/.*$/gm, '$1');
  const found: string[] = [];
  for (const m of text.matchAll(/\.subscribe\(/g)) {
    // The statement the call is in: back to the previous `;`, `{` or `}`.
    const before = text.slice(0, m.index);
    const start = Math.max(before.lastIndexOf(';'), before.lastIndexOf('{'), before.lastIndexOf('}')) + 1;
    const statement = before.slice(start);
    const name = /^\s*const (\w+)\s*=/.exec(statement)?.[1];
    if (!name || !text.includes(`${name}.unsubscribe()`)) {
      found.push(`${before.split('\n').length}: ${statement.trim().split('\n')[0]}`);
    }
  }
  return found;
}

describe('subscriptions in components', () => {
  it('finds an unowned one and accepts a kept one', () => {
    expect(unowned('x.subscribe(f);')).toHaveLength(1);
    expect(unowned('const s = x.subscribe(f);')).toHaveLength(1);
    expect(unowned('a();\n// Why; a note.\nconst s = x.subscribe(f);\nonCleanup(() => s.unsubscribe());')).toEqual([]);
    expect(unowned('const s = x\n  .pipe(take(1))\n  .subscribe(f);\nonCleanup(() => s.unsubscribe());')).toEqual([]);
  });

  it('are bound, or kept and unsubscribed', () => {
    const offenders: string[] = [];
    for (const file of components(src)) {
      for (const at of unowned(fs.readFileSync(file, 'utf8'))) offenders.push(`${path.relative(src, file)}:${at}`);
    }
    expect(offenders).toEqual([]);
  });
});
