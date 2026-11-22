/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The render configuration's table (_spec/README.md "Render", Configuration) against `DefaultRenderConfig`: every item
// of the config has a row (a row may name a group, `heights`, or several items, `labelWidth`, `labelGap`; `.hitPx`
// continues the item before it), and every row names an item the config has. Why a test: the table was kept by hand
// and a slice's close found seventeen rows missing.

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { DefaultRenderConfig } from '../render/config.js';
import { DefaultResolveConfig } from '../resolve/model.js';

/** The items a Markdown table names in its first column, after the heading matching `heading`. */
function tableItems(md: string, heading: RegExp): string[] {
  const lines = md.split('\n');
  let i = lines.findIndex((l) => heading.test(l));
  if (i < 0) throw new Error(`no heading ${heading}`);
  while (i < lines.length && !lines[i].startsWith('|')) i++;
  const items: string[] = [];
  for (i += 2; i < lines.length && lines[i].startsWith('|'); i++) {
    const first = lines[i].split('|')[1];
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

/** The config's items: paths to its values that are no plain non-empty object. */
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

/** The config's items no row covers, and the rows that name no item. */
function tableDrift(items: readonly string[], rows: readonly string[]) {
  const covered = (item: string) => rows.some((r) => item === r || item.startsWith(`${r}.`));
  const named = (row: string) => items.some((i) => i === row || i.startsWith(`${row}.`));
  return { missing: items.filter((i) => !covered(i)), unknown: rows.filter((r) => !named(r)) };
}

describe('the Configuration table of RenderConfig', () => {
  it('has a row for every item of DefaultRenderConfig, and names none it lacks', () => {
    const md = readFileSync(fileURLToPath(new URL('../../_spec/README.md', import.meta.url)), 'utf8');
    const rows = tableItems(md.slice(md.indexOf('\n## Render\n')), /^### Configuration$/);
    expect(rows.length).toBeGreaterThan(40);
    // `values.look` is internal until the Values looks story's pick is built as the default (then it goes, here too)
    expect(
      tableDrift(
        configItems(DefaultRenderConfig).filter((i) => i !== 'values.look'),
        rows,
      ),
    ).toEqual({ missing: [], unknown: [] });
  });
});

describe('the Configuration table of ResolveConfig', () => {
  it('has a row for every item of DefaultResolveConfig, and names none it lacks', () => {
    const md = readFileSync(fileURLToPath(new URL('../../_spec/README.md', import.meta.url)), 'utf8');
    const rows = tableItems(md.slice(md.indexOf('\n## Resolve\n')), /^### Configuration$/);
    expect(tableDrift(configItems(DefaultResolveConfig), rows)).toEqual({ missing: [], unknown: [] });
  });
});
