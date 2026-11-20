/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's Configuration table (_spec/README.md "Configuration") against `DefaultViewerConfig`: every item has a
// row and every row names an item (rows marked "in `SeqViewerSpec`" are the composition's parts, not config). Why a
// test: the table was kept by hand and a slice's close found rows missing.

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { DefaultViewerConfig } from '../config.js';

function tableItems(md: string, heading: RegExp): string[] {
  const lines = md.split('\n');
  let i = lines.findIndex((l) => heading.test(l));
  if (i < 0) throw new Error(`no heading ${heading}`);
  while (i < lines.length && !lines[i].startsWith('|')) i++;
  const items: string[] = [];
  for (i += 2; i < lines.length && lines[i].startsWith('|'); i++) {
    const first = lines[i].split('|')[1];
    if (first.includes('SeqViewerSpec')) continue;
    let parent = '';
    for (const m of first.matchAll(/`([^`]+)`/g)) {
      let item = m[1];
      if (item.startsWith('.')) item = parent + item;
      else parent = item.includes('.') ? item.slice(0, item.lastIndexOf('.')) : '';
      items.push(item);
    }
  }
  return items;
}

function configItems(config: object, prefix = ''): string[] {
  const out: string[] = [];
  for (const [k, v] of Object.entries(config)) {
    const path = prefix + k;
    if (v && typeof v === 'object' && !Array.isArray(v) && Object.keys(v).length)
      out.push(...configItems(v, `${path}.`));
    else out.push(path);
  }
  return out;
}

describe('the Configuration table of ViewerConfig', () => {
  it('has a row for every item of DefaultViewerConfig, and names none it lacks', () => {
    const md = readFileSync(fileURLToPath(new URL('../../_spec/README.md', import.meta.url)), 'utf8');
    const rows = tableItems(md, /^## Configuration$/);
    const items = configItems(DefaultViewerConfig);
    expect(rows.length).toBeGreaterThan(60);
    const covered = (item: string) => rows.some((r) => item === r || item.startsWith(`${r}.`));
    const named = (row: string) => items.some((i) => i === row || i.startsWith(`${row}.`));
    expect({ missing: items.filter((i) => !covered(i)), unknown: rows.filter((r) => !named(r)) }).toEqual({
      missing: [],
      unknown: [],
    });
  });
});
