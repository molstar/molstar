/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// KS1: keyed sets in the viewer (seqspec:format#layer-types): a rows layer ordered and filtered by key tables, computed
// ones among them (`row_stats`, `redundancy`), through the graph: the lanes are the resolved permutation, the layout's
// display rows follow them (`dataRow`), a recomputed table reorders them, the tree is off while the order is another's.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { set } from './axis-fixtures.js';
import { opened } from './layout-fixtures.js';

const ROWS = ['ACDEFGHIKL', 'ACDEFGHIKM', 'ACDEFGHIAM', 'ACDEFGHIKV', 'WYWYWYWYWY', 'WYWYWY----'];

function doc(rows: Record<string, unknown>) {
  const b = createSeqViewSpec({ title: 'Keyed' });
  const aln = b.data.alignment({ ref: 'aln', value: ROWS, labels: ROWS.map((_, i) => `r${i}`) });
  aln.tree({ value: '(((r0,r1),(r2,r3)),(r4,r5));' } as never);
  b.data.computed({ ref: 'stats', provider: 'row_stats', alignment: aln } as never);
  b.data.computed({
    ref: 'nr',
    provider: 'redundancy',
    alignment: aln,
    threshold: 0.8,
    order: { data: 'stats', field: 'identity', direction: 'desc' },
  } as never);
  const v = b.view({ ref: 'v', label: 'V', style: 'msa' }).block({ alignment: aln });
  v.track({ label: 'Rows' }).rows({ ref: 'members', ...rows } as never);
  return b.getDocument();
}
const rowsOf = (t: Awaited<ReturnType<typeof opened>>) => {
  const r = t.ctx.layout.visibleRows('v.main')?.current;
  return r ? Array.from(r.dataRow.subarray(0, r.count)) : [];
};

describe('a rows layer ordered by key tables', () => {
  it("orders by identity (desc), clusters, and keeps representatives only; each is the layout's display rows", async () => {
    const byId = await opened(doc({ order: { data: 'stats', field: 'identity', direction: 'desc' } }));
    expect(rowsOf(byId)).toEqual([0, 1, 3, 2, 4, 5]);
    const nr = await opened(
      doc({
        order: [
          { data: 'nr', field: 'cluster' },
          { data: 'stats', field: 'identity', direction: 'desc' },
        ],
      }),
    );
    expect(rowsOf(nr)).toEqual([0, 1, 3, 2, 4, 5]);
    const rep = await opened(doc({ filter: [{ data: 'nr', field: 'representative', values: [true] }] }));
    expect(rowsOf(rep)).toEqual([0, 4]);
    expect(rep.layout('v.main').rows.count).toBe(2);
  });

  it('a changed key table reorders the lanes: an edit to the alignment recomputes the statistics', async () => {
    const t = await opened(doc({ order: { data: 'stats', field: 'identity', direction: 'desc' } }));
    expect(rowsOf(t)).toEqual([0, 1, 3, 2, 4, 5]);
    // Row 4 becomes the nearest to the reference, row 1 the farthest of the first four.
    const next = ['ACDEFGHIKL', 'WYDEFGHIKM', 'ACDEFGHIAM', 'ACDEFGHIKV', 'ACDEFGHIKL', 'WYWYWY----'];
    await t.ctx.applyPatch(set('aln', { value: next }));
    expect(rowsOf(t)).toEqual([0, 4, 3, 2, 1, 5]);
  });

  it("a patch to the layer's order and filter reorders the rows (the story's buttons), and back", async () => {
    const t = await opened(doc({ order: 'data' }));
    expect(rowsOf(t)).toEqual([0, 1, 2, 3, 4, 5]);
    const byId = { data: 'stats', field: 'identity', direction: 'desc' };
    await t.ctx.applyPatch(
      set('members', { order: byId, filter: [{ data: 'nr', field: 'representative', values: [true] }] }),
    );
    expect(rowsOf(t)).toEqual([0, 4]);
    await t.ctx.applyPatch(set('members', { order: [{ data: 'nr', field: 'cluster' }, byId], filter: [] }));
    expect(rowsOf(t)).toEqual([0, 1, 3, 2, 4, 5]);
    await t.ctx.applyPatch(set('members', { order: 'data', filter: [] }));
    expect(rowsOf(t)).toEqual([0, 1, 2, 3, 4, 5]);
  });

  it("the tree is drawn while the order is the tree's, and folds nothing under another order", async () => {
    const tree = await opened(doc({ order: 'tree' }));
    expect(rowsOf(tree)).toEqual([0, 1, 2, 3, 4, 5]);
    const clade = 1;
    tree.ctx.ui.toggleClade('aln', clade);
    expect(tree.layout('v.main').rows.node).not.toBeNull();
    const keyed = await opened(doc({ order: { data: 'stats', field: 'identity', direction: 'desc' } }));
    keyed.ctx.ui.toggleClade('aln', clade);
    expect(keyed.layout('v.main').rows.node).toBeNull();
    expect(rowsOf(keyed)).toEqual([0, 1, 3, 2, 4, 5]);
  });
});
