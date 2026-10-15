/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Documents and loci for the composer's overlay (CS1b), with no Node import: the browser suite reads them too.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { type Loci, OrderedSet, type SeqSet } from '@molstar/seq-core/model';

/** Two proteins in one section with letters and values; stacked, wrapped (`text`) or side by side. */
export function twoBlocks(o: { stacked?: boolean; wrapped?: boolean } = {}): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEK'.repeat(12) });
  const c = b.data.sequence({ ref: 'c', value: 'GSWQERTYIP'.repeat(9) });
  const v = b.view({ ref: 'v', label: 'V', height: 'box', ...(o.wrapped ? { layout: 'text' } : {}) });
  v.block({ sequence: a });
  v.block({ sequence: c });
  v.track({ label: 'L' }).letters();
  const doc = b.getDocument();
  type Node = { kind: string; params?: Record<string, unknown>; children?: Node[] };
  const walk = (n: Node): void => {
    if (n.kind === 'section' && o.stacked) n.params = { ...n.params, arrangement: 'stacked' };
    n.children?.forEach(walk);
  };
  walk(doc.root as Node);
  return doc;
}

/** An alignment's rows ("MK-VLA": a gap at column 2), a data column, and a pinned ruler unless `bare`. */
export function family(o: { bare?: boolean } = {}): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const sto = ['# STOCKHOLM 1.0', 'a1 MK-VLA', 'a2 MKVQLA', 'a3 MRV-LA', 'a4 MKVQLA', '//', ''].join('\n');
  const aln = b.data.alignment({ ref: 'aln', source: { value: sto, format: 'stockholm' } });
  const meta = b.data.annotation({
    ref: 'meta',
    geometry: 'rows',
    target: aln,
    value: { row: [0, 1, 2, 3], score: [1, 2, 3, 4] },
  } as never);
  const v = b.view({ ref: 'f', label: 'F', height: 'box' });
  v.block({ alignment: aln });
  if (!o.bare) v.track({ label: '' }).ruler();
  v.track({ label: 'Rows' }).rows({ data_columns: [{ data: meta, field: 'score', label: 'Score' }] } as never);
  return b.getDocument();
}

/** Seq loci on a sequence's unit over [s, e). */
export function seqLoci(set: SeqSet, sequence: string, s: number, e: number): Loci {
  const unit = set.units.key.findIndex((k, u) => set.units.alive[u] === 1 && set.strings.get(k) === `${sequence}/unit`);
  return { kind: 'seq-loci', set, elements: [{ unit, indices: OrderedSet.ofRange(s, e) }] } as unknown as Loci;
}
/** Row loci on alignment 0's row over residues [s, e). */
export function rowLoci(set: SeqSet, row: number, s: number, e: number): Loci {
  const elements = [{ row, indices: OrderedSet.ofRange(s, e) }];
  return { kind: 'row-loci', set, alignment: 0, elements } as unknown as Loci;
}
