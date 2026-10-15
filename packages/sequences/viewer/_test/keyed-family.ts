/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The keyed-sets story's document and its five modes (`Seq*/Views and layouts/Alignments › Order and filter rows`), shared with the test that applies each mode's
// patch: PF00042's rows with `row_stats` and two clusterings as data columns, and the `rows` layer's `order` and `filter`
// the buttons set. The layer's ref is the patch's target (a track's ref is no layer's).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';

export interface FamilyRow {
  readonly label: string;
  readonly gapped: string;
}
/** The rows layer a mode's patch sets. */
export const MEMBERS_LAYER = 'members.rows';

const BY_IDENTITY = { data: 'stats', field: 'identity', direction: 'desc' };
export interface Mode {
  readonly key: string;
  readonly label: string;
  readonly order: unknown;
  readonly filter: readonly unknown[];
}
/** The layer's `order` and `filter` for each choice; the buttons patch the layer's params. */
export const MODES: readonly Mode[] = [
  { key: 'data', label: 'Data order', order: 'data', filter: [] },
  { key: 'identity', label: 'By identity to HBA_HUMAN', order: BY_IDENTITY, filter: [] },
  {
    key: 'clusters',
    label: 'Clustered at 90%, then by identity',
    order: [{ data: 'nr90', field: 'cluster' }, BY_IDENTITY],
    filter: [],
  },
  {
    key: 'reps90',
    label: 'Representatives at 90%',
    order: BY_IDENTITY,
    filter: [{ data: 'nr90', field: 'representative', values: [true] }],
  },
  {
    key: 'reps50',
    label: 'Representatives at 50%',
    order: BY_IDENTITY,
    filter: [{ data: 'nr50', field: 'representative', values: [true] }],
  },
];

/** The patch that sets a mode on the members layer. */
export const modePatch = (m: Mode) => ({
  format: 'seq_view_spec_patch' as const,
  ops: [{ op: 'set' as const, ref: MEMBERS_LAYER, params: { order: m.order, filter: m.filter } }],
});

/** The family (32 globins) with `row_stats` and two clusterings as data columns, in data order. */
export function keyedFamilyDoc(family: readonly FamilyRow[]): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Keyed sets on PF00042' });
  const aln = b.data.alignment({
    ref: 'pf00042',
    value: family.map((r) => r.gapped),
    labels: family.map((r) => r.label),
  });
  const stats = b.data.computed({ ref: 'stats', provider: 'row_stats', alignment: aln } as never);
  // Clusters are taken from the most to the least identical to the reference row, so the reference leads its own.
  const cluster = (ref: string, threshold: number) =>
    b.data.computed({ ref, provider: 'redundancy', alignment: aln, threshold, order: BY_IDENTITY } as never);
  const nr90 = cluster('nr90', 0.9);
  const nr50 = cluster('nr50', 0.5);
  const v = b
    .view({ ref: 'v', label: 'PF00042', style: 'msa', height: 'box', height_rows: 18 })
    .block({ alignment: aln });
  v.track({ label: '' }).ruler();
  v.track({ ref: 'members', label: 'Members' }).rows({
    ref: MEMBERS_LAYER,
    order: 'data',
    data_columns: [
      { data: stats, field: 'identity', label: 'Identity', type: 'fraction' },
      { data: stats, field: 'coverage', label: 'Coverage', type: 'fraction' },
      { data: nr90, field: 'cluster', label: 'Cluster 90%', type: 'number' },
      { data: nr50, field: 'cluster', label: 'Cluster 50%', type: 'number' },
    ],
  } as never);
  return b.getDocument();
}
