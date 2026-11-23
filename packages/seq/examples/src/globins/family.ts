/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { matchLeaves, readNewick } from '@molstar/seqx-core/data/alignment-readers';
import {
  createSeqSpec,
  type NavigatorOf,
  type NavigatorParams,
  overviewNavigator,
  theme,
} from '@molstar/seqx-core/lang';
import { buildClades, cladesAtDepth } from '@molstar/seqx-core/render';
import { balancedNewick, familyRegistry, GeneratedFormat, SeedFile, SeedTree } from './data.js';
import { familySources } from './sources.js';

// The Globins views (EC-38 to EC-44, EC-97, EC-100, EC-103, EC-106, EC-154): the Pfam seed with its conservation logo
// and its rows by tree, each row's species (UniProt's, a file joined by row key) on the left and its identity to the
// reference row on the right; its clades folded to summary rows; the ESPript-style figure of its rows in text; the
// 10,000-row variant (held, by tree) and a lazy one of 10,000,000 rows (data order, its rows made on demand).

export const title = 'Globins';

export const registry = familyRegistry;

/** The seed as published (read by location), its tree inline, the rows' columns file joined by row key. */
function familyData(b: ReturnType<typeof createSeqSpec>) {
  const { aln } = familySources(b);
  const tree = b.data.tree('tree', { of: aln, newick: SeedTree });
  return { aln, tree };
}

/** The family's row columns: the species on the left, the identity to the reference row on the right. */
export const FamilyColumns = [
  { field: 'OS', title: 'Species', width: 104 },
  { field: 'identity', side: 'right', as: 'number', title: 'Identity' },
] as const;

/**
 * The vertical overview of an alignment's rows (seqx-core _spec "Render", Navigators): the rows summarized down a
 * column at the right, its window the rows and the columns in view, a grip each side (top and bottom for the rows, left
 * and right for the columns), 96 px wide.
 */
export function rowsOverview(of: NavigatorOf): NavigatorParams {
  return { ...overviewNavigator(of), axes: 'both', size: OverviewWidth };
}

/** The vertical overview's width, px (seqx-core _spec "Render", Navigators). */
export const OverviewWidth = 96;

export interface FamilyOptions {
  readonly collapsed?: readonly string[];
  /** The vertical overview of the rows at the right (`true` when absent; without it the rows rail). */
  readonly overview?: boolean;
  /** The rows' residue scheme (the layer's default when absent). */
  readonly scheme?: string;
}

/** The conservation logo's track height, px: its letters read at the fit (the default is the thin strip a value track gets). */
export const LogoHeight = 64;

/** The family: a ruler, its conservation logo, its rows by tree beside the tree, species and identity beside them. */
export function build(options: FamilyOptions = {}): Node {
  const b = createSeqSpec({ title: 'Globins: the family' });
  const { aln, tree } = familyData(b);
  b.views.view('main', (v) => {
    const s = v.section('family', { of: aln });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    s.track('logo', { label: 'Conservation', height: LogoHeight }).layer('logo', { type: 'logo', data: aln });
    s.track('rows', { label: '' }).layer('rows', {
      type: 'rows',
      data: aln,
      tree,
      order: 'tree',
      columns: FamilyColumns.map((c) => ({ ...c })),
      ...(options.collapsed ? { collapsed: [...options.collapsed] } : {}),
      ...(options.scheme ? { color: theme.scheme(options.scheme) } : {}),
    });
    if (options.overview !== false) v.navigator('overview', rowsOverview(s));
  });
  return b.build();
}

/** The family with clades folded to summary rows (their consensus, "first +n"): the opening state a writer kept. */
export function folded(): Node {
  return build({ collapsed: cladesToFold() });
}

/**
 * Two clades of three rows or more at depth 2 of the seed's tree, by their names (`first..last`): what the folded view
 * opens with. The tree's own leaves stand for the rows (a clade's name is its leaves', whatever the rows' order).
 */
export function cladesToFold(): string[] {
  const tree = readNewick(SeedTree);
  const rows = tree.leaves.map((n) => ({ label: tree.label[n] }));
  const clades = buildClades(tree, matchLeaves(tree, rows), rows.length)!;
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

/** The scale variant held (10,000 × 300 by default), by its balanced tree; its rows sampled from the seed's profile. */
export function scale(
  options: { readonly rows?: number; readonly columns?: number; readonly overview?: boolean } = {},
): Node {
  const rows = options.rows ?? 10_000;
  const b = createSeqSpec({ title: `Globins: ${rows.toLocaleString('en-US')} rows` });
  const src = b.data.source('generated', {
    format: GeneratedFormat,
    value: { rows, columns: options.columns ?? 300, seed: 1, profile: SeedFile },
  });
  const aln = b.data.alignment('scale', { alphabet: 'protein', source: src });
  const tree = b.data.tree('tree', { of: aln, newick: balancedNewick(rows) });
  b.views.view('main', (v) => {
    const s = v.section('rows', { of: aln });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    s.track('logo', { label: 'Conservation', height: LogoHeight }).layer('logo', { type: 'logo', data: aln });
    s.track('rows', { label: '' }).layer('rows', { type: 'rows', data: aln, tree, order: 'tree' });
    if (options.overview !== false) v.navigator('overview', rowsOverview(s));
  });
  return b.build();
}

/** The lazy variant: 10,000,000 rows in data order, only the rows in view ever made (EC-41, EC-44). */
export function lazy(options: { readonly rows?: number } = {}): Node {
  const rows = options.rows ?? 10_000_000;
  const b = createSeqSpec({ title: `Globins: ${rows.toLocaleString('en-US')} rows, lazy` });
  const src = b.data.source('generated', {
    format: GeneratedFormat,
    value: { rows, columns: 300, seed: 2, profile: SeedFile },
  });
  const aln = b.data.alignment('lazy', { alphabet: 'protein', source: src });
  b.views.view('main', (v) => {
    const s = v.section('rows', { of: aln });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    s.track('rows', { label: '' }).layer('rows', { type: 'rows', data: aln });
  });
  return b.build();
}
