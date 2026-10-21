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
import { extensionPrismFiles } from './extension-views.js';

// Nothing in Prism observes DOM mutations over a large subtree (prism AGENTS.md: views render what their states say;
// a MutationObserver over a panel is told of every row the lists make and every class a render sets). The few that
// watch a small subtree of their own are listed here with why; any other `new MutationObserver(` in Prism's source (ui,
// viewer, assistant, the extensions' Prism UI; not tests, stories or testing aids) fails.

const Allowed: Record<string, string> = {
  // `rovingTabs`: a tab list's own tabs (a handful), to keep one Tab stop as tabs come, go and are selected.
  'packages/prism/ui/src/dom/focus.ts': 'a tab list, its tabs',
};

const prism = path.resolve(import.meta.dirname, '../../..');
const repo = path.resolve(prism, '../..');

function* sources(dir: string): Generator<string> {
  if (!fs.existsSync(dir)) return;
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (['_test', '_stories', '_generated', 'testing', 'node_modules', 'lib'].includes(entry.name)) continue;
      yield* sources(full);
    } else if (/\.tsx?$/.test(entry.name) && !entry.name.endsWith('.d.ts')) yield full;
  }
}

describe('no MutationObserver over the UI', () => {
  it('none in Prism but the listed ones', () => {
    const found: string[] = [];
    const scan = (file: string) => {
      const text = fs.readFileSync(file, 'utf8').replace(/(^|[^:])\/\/.*$/gm, '$1');
      if (/new\s+MutationObserver\s*\(/.test(text)) found.push(path.relative(repo, file).split(path.sep).join('/'));
    };
    for (const pkg of ['ui', 'viewer', 'assistant']) for (const f of sources(path.join(prism, pkg, 'src'))) scan(f);
    for (const f of extensionPrismFiles('all')) scan(f);
    expect(found.filter((f) => !(f in Allowed))).toEqual([]);
    // A listed one gone: take it off the list.
    expect(Object.keys(Allowed).filter((f) => !found.includes(f))).toEqual([]);
  });
});
