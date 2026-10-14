/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A synthetic protein family for vertical scrolling, the overview and row zoom (synthetic: not real sequences): the
// SH2 excerpt's rows (case 3) as seeds, each member a seed with deterministic substitutions and an occasional gap.
// No file, no network: an alignment by value, a box section of `height_rows` rows under a ruler and a conservation track.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';

const SEEDS = [
  'WYFGKITRRESERLLLNAENPRGTFLVRES',
  'WFFGKIPRAKAEEMLSKQRHDGAFLIRESE',
  'WFHPNITGVEKDAENLLLTRGVDGSFLARP',
  'WYWGDISREEVNEKLRDTADGTFLVRDASS',
];
const AMINO = 'ACDEFGHIKLMNPQRSTVWY';

/** The family's rows: `n` members over the seeds' 30 columns. */
export function familyRows(n: number): { labels: string[]; rows: string[] } {
  const labels: string[] = [];
  const rows: string[] = [];
  for (let k = 0; k < n; k++) {
    const seed = SEEDS[k % SEEDS.length];
    const row = [...seed].map((c, i) => {
      const h = (k * 131 + i * 17 + ((k * i) % 7)) % 100;
      return h < 12 ? AMINO[(k + i * 3) % AMINO.length] : h < 15 ? '-' : c;
    });
    labels.push(`SYN${String(k + 1).padStart(3, '0')}/1-30`);
    rows.push(row.join(''));
  }
  return { labels, rows };
}

export const FAMILY_SECTION = 'fam.main';

/** A box section of `heightRows` rows over an `n`-row synthetic family, in the `msa` preset. */
export function syntheticFamily(n = 120, heightRows = 12): SeqViewSpecDoc {
  const { labels, rows } = familyRows(n);
  const width = Math.max(...labels.map((l) => l.length)) + 2;
  const sto = ['# STOCKHOLM 1.0', '#=GF ID   SYN-SH2', ...rows.map((r, k) => labels[k].padEnd(width) + r), '//', ''];
  const b = createSeqViewSpec({ title: `A synthetic SH2-like family of ${n} rows` });
  const aln = b.data.alignment({ ref: 'family', source: { value: sto.join('\n'), format: 'stockholm' } });
  const cons = b.data.computed({ ref: 'cons', provider: 'column_stats', alignment: aln, measure: 'conservation' });
  const v = b.view({ ref: 'fam', label: 'Family', style: 'msa', height: 'box', height_rows: heightRows });
  v.block({ alignment: aln });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Conservation' }).values({ data: cons });
  v.track({ label: 'Members' }).rows();
  return b.getDocument();
}

/**
 * The synthetic family with a synthetic tree (deterministic, not a real phylogeny): a clade per seed (members k ≡ s
 * mod 4), each a balanced binary tree over its members with branch lengths by depth, in Newick by value. Optional data
 * columns from a `rows` annotation (a number, a fraction, a text, a category per member) and a points layer on rows
 * (two sites, an active and a binding site, so its legend names them).
 */
export function familyWithTree(n: number, heightRows: number, o: { columns?: boolean; sites?: boolean } = {}) {
  const { labels, rows } = familyRows(n);
  const clade = (ks: number[], depth: number): string =>
    ks.length === 1
      ? `'${labels[ks[0]]}':${(0.05 + (ks[0] % 3) * 0.03).toFixed(2)}`
      : `(${clade(ks.slice(0, ks.length >> 1), depth + 1)},${clade(ks.slice(ks.length >> 1), depth + 1)}):${(0.3 / (depth + 1)).toFixed(2)}`;
  const seeds = [0, 1, 2, 3].map((s) => Array.from({ length: n }, (_, k) => k).filter((k) => k % 4 === s));
  const newick = `((${clade(seeds[0], 1)},${clade(seeds[1], 1)}):0.4,(${clade(seeds[2], 1)},${clade(seeds[3], 1)}):0.4);`;
  const width = Math.max(...labels.map((l) => l.length)) + 2;
  const sto = ['# STOCKHOLM 1.0', ...rows.map((r, k) => labels[k].padEnd(width) + r), '//', ''];
  const b = createSeqViewSpec({ title: `A synthetic family of ${n} rows with a synthetic tree` });
  const aln = b.data.alignment({ ref: 'family', source: { value: sto.join('\n'), format: 'stockholm' } });
  aln.tree({ source: { value: newick, format: 'newick' } });
  const v = b.view({ ref: 'fam', label: 'Family', style: 'msa', height: 'box', height_rows: heightRows });
  v.block({ alignment: aln });
  v.track({ label: '' }).ruler();
  const members = v.track({ label: 'Members' });
  const meta = o.columns
    ? b.data.annotation({
        ref: 'meta',
        geometry: 'rows',
        target: aln,
        value: {
          row: Array.from({ length: n }, (_, k) => k),
          length: Array.from({ length: n }, (_, k) => 28 + (k % 3)),
          identity: Array.from({ length: n }, (_, k) => ((k * 37) % 100) / 100),
          species: Array.from({ length: n }, (_, k) => `Species ${k % 7}`),
          clade: Array.from({ length: n }, (_, k) => `clade ${k % 4}`),
        },
      })
    : null;
  members.rows(
    meta
      ? {
          data_columns: [
            { data: meta, field: 'length', label: 'Len', type: 'number' },
            { data: meta, field: 'identity', label: 'Id.', type: 'fraction' },
            { data: meta, field: 'species', label: 'Species', type: 'text' },
            { data: meta, field: 'clade', label: 'Clade', type: 'category' },
          ],
        }
      : {},
  );
  if (o.sites) {
    const sites = b.data.annotation({
      ref: 'sites',
      geometry: 'ranges',
      target: aln,
      numbering: 'ordinal',
      value: { target: [aln.row(0), aln.row(5)], start: [3, 10], category: ['Active site', 'Binding site'] },
    } as never);
    members.points({ data: sites });
  }
  return b.getDocument();
}

/** A family `n` rows by 30 × `copies` columns (the fit cell is small enough to reach every level of detail). */
export function wideFamily(n: number, copies: number) {
  const { labels, rows } = familyRows(n);
  const sto = ['# STOCKHOLM 1.0', ...rows.map((r, k) => `${labels[k].padEnd(16)}${r.repeat(copies)}`), '//', ''];
  const b = createSeqViewSpec({ title: 'A wide family' });
  const aln = b.data.alignment({ ref: 'family', source: { value: sto.join('\n'), format: 'stockholm' } });
  const v = b.view({ ref: 'fam', label: 'Family', style: 'msa', height: 'box', height_rows: 12 });
  v.block({ alignment: aln });
  v.track({ label: 'Members' }).rows();
  return b.getDocument();
}
