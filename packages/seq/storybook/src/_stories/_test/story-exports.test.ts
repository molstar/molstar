/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Every named export of a stories file is a story: Storybook renders each one, so an exported table (a map of where
// the views are) is a broken story ("the component annotation is missing"). Read as text, so no story is imported.

import { readdirSync, readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');

function storiesFiles(dir: string): string[] {
  const out: string[] = [];
  for (const e of readdirSync(dir, { withFileTypes: true })) {
    if (e.name === 'node_modules' || e.name === '_test') continue;
    const p = join(dir, e.name);
    if (e.isDirectory()) out.push(...storiesFiles(p));
    else if (e.name.endsWith('.stories.ts')) out.push(p);
  }
  return out;
}

describe('stories files', () => {
  it('export stories only: an annotated export is a StoryObj', () => {
    const files = storiesFiles(root);
    expect(files.length).toBeGreaterThan(20);
    const wrong: string[] = [];
    for (const f of files) {
      for (const m of readFileSync(f, 'utf8').matchAll(/^export const (\w+)\s*:\s*([^=]+)=/gm)) {
        if (!m[2].trim().startsWith('StoryObj')) wrong.push(`${f.slice(root.length + 1)}: ${m[1]}`);
      }
    }
    expect(wrong).toEqual([]);
  });

  it('export no function or class: Storybook lists it as a story too', () => {
    const wrong: string[] = [];
    for (const f of storiesFiles(root)) {
      for (const m of readFileSync(f, 'utf8').matchAll(/^export (?:async )?(?:function\*?|class) (\w+)/gm))
        wrong.push(`${f.slice(root.length + 1)}: ${m[1]}`);
    }
    expect(wrong).toEqual([]);
  });
});
