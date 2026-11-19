/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { matchLeaves, readNewick } from '@molstar/seqx-core/data/alignment-readers';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { buildClades, cladesAtDepth } from '@molstar/seqx-core/render';
import { balancedNewick, GeneratedFormat, globinFamily } from './data.js';

// The Globins views (EC-38 to EC-44, EC-97, EC-100, EC-103, EC-106, EC-154): the family with its conservation logo and
// its rows by tree; its clades folded to summary rows; the ESPript-style figure of its rows in text; the 10,000-row
// variant (held, by tree) and a lazy one of 10,000,000 rows (data order, its rows made on demand).

export const title = 'Globins';

/** The family's alignment and tree, inline. */
function familyData(b: ReturnType<typeof createSeqSpec>) {
  const f = globinFamily();
  const aln = b.data.alignment('globins', {
    alphabet: 'protein',
    rows: { label: [...f.label], letters: [...f.letters] },
    reference: f.label[0],
  });
  const tree = b.data.tree('tree', { of: aln, newick: f.newick });
  return { aln, tree, family: f };
}

/** The family: a ruler, its conservation logo, its rows by tree beside the tree. */
export function build(options: { readonly collapsed?: readonly string[] } = {}): Node {
  const b = createSeqSpec({ title: 'Globins: the family' });
  const { aln, tree } = familyData(b);
  b.views.view('main', (v) => {
    const s = v.section('family', { of: aln });
    s.track('ruler', { label: 'Columns' }).layer('ruler', { type: 'ruler' });
    s.track('logo', { label: 'Conservation' }).layer('logo', { type: 'logo', data: aln });
    s.track('rows', { label: '' }).layer('rows', {
      type: 'rows',
      data: aln,
      tree,
      order: 'tree',
      ...(options.collapsed ? { collapsed: [...options.collapsed] } : {}),
    });
  });
  return b.build();
}

/** The family with clades folded to summary rows (their consensus, "first +n"): the opening state a writer kept. */
export function folded(): Node {
  const f = globinFamily();
  // the tree's two largest clades under the root's first child, by name (`first..last`)
  return build({ collapsed: cladesToFold(f) });
}

/** Two clades of three rows or more at depth 2 of the family's tree, by their names (what the folded view opens with). */
function cladesToFold(f: ReturnType<typeof globinFamily>): string[] {
  const tree = readNewick(f.newick);
  const clades = buildClades(
    tree,
    matchLeaves(
      tree,
      f.label.map((label) => ({ label })),
    ),
    f.label.length,
  )!;
  return cladesAtDepth(clades, 2)
    .filter((c) => clades.hi[c] - clades.lo[c] >= 3)
    .slice(0, 2)
    .map((c) => clades.name[c]);
}

/** The ESPript-style figure: the family's rows as text, a ruler over them, every row on every line. */
export function figure(): Node {
  const b = createSeqSpec({ title: 'Globins: the alignment figure' });
  const { aln } = familyData(b);
  b.views.view('main', (v) => {
    const s = v.section('figure', { of: aln, layout: 'text', group: 10 });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    s.track('rows', { label: '' }).layer('rows', { type: 'rows', data: aln });
  });
  return b.build();
}

/** The scale variant held (10,000 × 300 by default), by its balanced tree. */
export function scale(options: { readonly rows?: number; readonly columns?: number } = {}): Node {
  const rows = options.rows ?? 10_000;
  const b = createSeqSpec({ title: `Globins: ${rows.toLocaleString('en-US')} rows` });
  const src = b.data.source('generated', {
    format: GeneratedFormat,
    value: { rows, columns: options.columns ?? 300, seed: 1 },
  });
  const aln = b.data.alignment('scale', { alphabet: 'protein', source: src });
  const tree = b.data.tree('tree', { of: aln, newick: balancedNewick(rows, options.columns ?? 300, 1) });
  b.views.view('main', (v) => {
    const s = v.section('rows', { of: aln });
    s.track('ruler', { label: 'Columns' }).layer('ruler', { type: 'ruler' });
    s.track('logo', { label: 'Conservation' }).layer('logo', { type: 'logo', data: aln });
    s.track('rows', { label: '' }).layer('rows', { type: 'rows', data: aln, tree, order: 'tree' });
  });
  return b.build();
}

/** The lazy variant: 10,000,000 rows in data order, only the rows in view ever made (EC-41, EC-44). */
export function lazy(options: { readonly rows?: number } = {}): Node {
  const rows = options.rows ?? 10_000_000;
  const b = createSeqSpec({ title: `Globins: ${rows.toLocaleString('en-US')} rows, lazy` });
  const src = b.data.source('generated', { format: GeneratedFormat, value: { rows, columns: 300, seed: 2 } });
  const aln = b.data.alignment('lazy', { alphabet: 'protein', source: src });
  b.views.view('main', (v) => {
    const s = v.section('rows', { of: aln });
    s.track('ruler', { label: 'Columns' }).layer('ruler', { type: 'ruler' });
    s.track('rows', { label: '' }).layer('rows', { type: 'rows', data: aln });
  });
  return b.build();
}
