/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Selection across changes (seqspec:interaction#selection; ported from the managers test) and held loci carried within
// a lineage (VS3); refs to loci with range args (seqspec:interaction#focus), item ids as a set.

import {
  type ItemLoci,
  type MoleculeInput,
  OrderedSet,
  type SeqLoci,
  type SeqSet,
  SeqSets,
} from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { createSeqViewer, type SeqViewerContext, type WritableSeqState } from '../src/index.js';
import { DefaultSeqViewerSpec } from '../src/default-spec.js';
import { carry } from '../src/lineage.js';
import { patch, smallDoc, viewer } from './graph-fixtures.js';

const unitKeys = (ctx: SeqViewerContext) => {
  const l = ctx.selection.current.value.seq as SeqLoci;
  return l.elements.map((e) => l.set.strings.get(l.set.units.key[e.unit]));
};
const mol = (key: string, residues: string): MoleculeInput => ({ key, residues, ownUnit: true });
const p53 = mol('p53', 'MEEPQSDPSVEPPLSQETFSDLWKLL');
const q = mol('q', 'MKVLATGHEK');
const hot = (stop: number) => ({
  key: 'hot',
  geometry: 'ranges' as const,
  columns: { target: ['p53/unit'], start: [3], stop: [stop], id: ['h1'] },
});
/** An item loci's extents through its set's tables: [start, stop) per row of its items. */
const extents = (l: ItemLoci) => {
  const ann = l.set.annotations[l.elements[0].annotation];
  const out: [number, number][] = [];
  for (let r = 0; ann && r < ann.count; r++)
    if (OrderedSet.has(l.elements[0].items, ann.item[r])) out.push([ann.start[r], ann.stop[r]]);
  return out;
};

describe('selection across changes', () => {
  it('kept as is across a patch that keeps the set', async () => {
    const { ctx } = viewer();
    await ctx.setDocument(smallDoc());
    ctx.selection.set(ctx.refs.toLoci('a', { start: 2, stop: 5 }), 'api');
    const before = ctx.selection.current.value;
    await ctx.applyPatch(patch([{ op: 'set', ref: 'v.main', params: { positions: 'dom' } }]));
    expect(ctx.selection.current.value).toBe(before);
  });

  it('carried across setData with build inputs; elements on tombstoned ids dropped', async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    await ctx.setData({ inputs: { molecules: [p53, q], annotations: [hot(6)] } });
    ctx.selection.set(ctx.refs.toLoci('p53/unit', { start: 0, stop: 3 }), 'api');
    ctx.selection.add(ctx.refs.toLoci('q/unit', { start: 1, stop: 2 }) as SeqLoci, 'pointer');
    const before = ctx.selection.current.value.seq as SeqLoci;
    expect(unitKeys(ctx)).toEqual(['p53/unit', 'q/unit']);
    const lineage = ctx.graph.assembled.value?.set.lineage;
    await ctx.setData({ inputs: { molecules: [p53, q], annotations: [hot(8)] } });
    // The same owners: rebound to the new set with the same indices.
    const kept = ctx.selection.current.value.seq as SeqLoci;
    expect(kept.set).toBe(ctx.graph.assembled.value?.set);
    expect(kept.elements).toEqual(before.elements);
    await ctx.setData({ inputs: { molecules: [p53], annotations: [hot(8)] } });
    expect(ctx.graph.assembled.value?.set.lineage).toBe(lineage);
    expect(unitKeys(ctx)).toEqual(['p53/unit']);
    expect(ctx.selection.current.value.origin).toBe('graph');
  });

  it('remapped by key across lineages, dropping what is gone', async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    const build = (ms: MoleculeInput[]) => SeqSets.build({ molecules: ms }).set;
    await ctx.setData({ set: build([p53, q]) });
    ctx.selection.set(ctx.refs.toLoci('q/unit', { start: 2, stop: 4 }), 'api');
    ctx.selection.add(ctx.refs.toLoci('p53/unit', { start: 0, stop: 1 }) as SeqLoci, 'api');
    const r = mol('r', 'MKV');
    await ctx.setData({ set: build([r, q]) }); // a new lineage: q's unit has another id, p53 is gone
    const sel = ctx.selection.current.value.seq as SeqLoci;
    expect(sel.set).toBe(ctx.graph.assembled.value?.set);
    expect(unitKeys(ctx)).toEqual(['q/unit']);
    expect(sel.elements[0].unit).toBe(1);
    expect(ctx.selection.current.value.origin).toBe('graph');
    expect(ctx.selection.isWhole(ctx.refs.toLoci('q/unit', { start: 2, stop: 4 }) as SeqLoci)).toBe(true);
  });
});

describe('held loci follow their inputs (VS3)', () => {
  it('an item selection after a range edit marks the new extent', async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    await ctx.setData({ inputs: { molecules: [p53], annotations: [hot(6)] } });
    ctx.selection.set(ctx.refs.toLoci('hot'), 'api');
    expect(extents(ctx.selection.current.value.items as ItemLoci)).toEqual([[3, 6]]);
    await ctx.setData({ inputs: { molecules: [p53], annotations: [hot(12)] } });
    const items = ctx.selection.current.value.items as ItemLoci;
    expect(items.set).toBe(ctx.graph.assembled.value?.set);
    expect(extents(items)).toEqual([[3, 12]]);
  });

  it('an assembled set that goes null and comes back keeps the selection, carried from the last set', async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    const build = (ms: MoleculeInput[]) => SeqSets.build({ molecules: ms }).set;
    await ctx.setData({ set: build([p53, q]) });
    ctx.selection.set(ctx.refs.toLoci('q/unit', { start: 2, stop: 4 }), 'api');
    const before = ctx.selection.current.value;
    (ctx.graph.assembled as WritableSeqState<unknown>).set(null);
    expect(ctx.selection.current.value).toBe(before);
    await ctx.setData({ set: build([mol('r', 'MKV'), q]) });
    expect(unitKeys(ctx)).toEqual(['q/unit']);
    expect((ctx.selection.current.value.seq as SeqLoci).set).toBe(ctx.graph.assembled.value?.set);
  });

  it('carry: the same set is the loci itself; nothing left is null', () => {
    const s1 = SeqSets.build({ molecules: [p53] }).set;
    const l: SeqLoci = {
      kind: 'seq-loci',
      set: s1,
      elements: [{ unit: 0 as never, indices: OrderedSet.ofRange(0, 3) as never }],
    };
    expect(carry(l, s1)).toBe(l);
    expect(carry(l, SeqSets.build({ molecules: [q] }).set as SeqSet)).toBe(null);
  });
});

describe('refs with item ids', () => {
  it('unknown, empty and trailing ids are dropped; none left is null', async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    const feats = {
      key: 'feats',
      geometry: 'ranges' as const,
      columns: {
        target: ['p53/unit', 'p53/unit', 'p53/unit'],
        start: [0, 4, 8],
        stop: [2, 6, 10],
        id: ['a', null, 'c'],
      },
    };
    await ctx.setData({ inputs: { molecules: [p53], annotations: [feats] } });
    const items = (args: readonly string[]) => {
      const l = ctx.refs.toLoci('feats', { items: args }) as ItemLoci | null;
      if (!l) return null;
      const out: number[] = [];
      OrderedSet.forEach(l.elements[0].items, (i) => void out.push(i), null);
      return out;
    };
    expect(items(['TYPO'])).toBe(null);
    expect(items([''])).toBe(null);
    expect(items(['a', ''])).toEqual([0]);
    expect(items(['c', 'a', 'a'])).toEqual([0, 2]);
  });
});
