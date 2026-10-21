/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';

// The extensions' Prism code, for the checks that hold for Prism's own (no rxjs in views, scoped attachments, owned
// subscriptions, no MutationObserver): an extension's views follow the same rules as the viewer's (prism AGENTS.md).
// Read by path, as data: prism-ui never imports an extension (opt-in-boundary.test.ts).

const Skipped = ['_test', '_stories', '_generated', 'testing', 'node_modules', 'lib'];

/** The repository's root */
export const RepoRoot = path.resolve(import.meta.dirname, '../../../../..');

function* sourceFiles(dir: string): Generator<string> {
  if (!fs.existsSync(dir)) return;
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (!Skipped.includes(entry.name)) yield* sourceFiles(full);
    } else if (/\.tsx?$/.test(entry.name) && !entry.name.endsWith('.d.ts')) yield full;
  }
}

/**
 * The source files of the extensions that use prism-ui (its name in the file: an import, the JSX runtime pragma),
 * `.tsx` only (the views) or all of them
 */
export function extensionPrismFiles(which: 'views' | 'all' = 'views'): string[] {
  const out: string[] = [];
  const extensions = path.join(RepoRoot, 'extensions');
  for (const ext of fs.readdirSync(extensions)) {
    for (const f of sourceFiles(path.join(extensions, ext, 'src'))) {
      if (which === 'views' && !f.endsWith('.tsx')) continue;
      if (fs.readFileSync(f, 'utf8').includes('@molstar/prism-ui')) out.push(f);
    }
  }
  return out;
}
