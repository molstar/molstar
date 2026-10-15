/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// NP1b: nodes a provider computes from the set (seqspec:model#alignments-and-trees, seqspec:format#mappings) in the
// transform graph: a `realign` alignment is assembled with the mappings, after the base set, so it reads its inputs
// there; a `codon_rows` mapping reads its two alignments. The graph's set is compile's, an edit to an input recomputes
// the alignment, and one to something else does not.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import {
  compile,
  createRegistries,
  KernelRegistry,
  MvsRegistry,
  register,
  type SeqViewSpecDoc,
} from '@molstar/seq-core/format';
import { type Alignment, KeyKind, type KeyTarget, type SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { Q } from '../../core/src/query/index.js';
import type { SeqViewerContext } from '../src/index.js';
import { patch, refsRan, viewer } from './graph-fixtures.js';

const S = ['MKVLAGWTQRSTEDNLK', 'MKVLAGWTRRSTEDNLK', 'KVLAGWTQRSTEDNLK'];
const CDNA = ['ATGGCTTCC', 'ATGAAAGTT'];

function realignDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  for (const [k, v] of S.entries()) b.data.sequence({ ref: `s${k}`, value: v });
  b.data.alignment({ ref: 're', provider: 'realign', inputs: ['s0', 's1', 's2'] } as never);
  b.data.alignment({ ref: 'prot', value: ['MAS', 'MKV'], alphabet: 'protein' });
  b.data.alignment({ ref: 'dna', value: CDNA, alphabet: 'dna' });
  b.data.mapping({ ref: 'cr', from: 'prot', to: 'dna', method: 'provider', provider: 'codon_rows' } as never);
  b.data.annotation({
    ref: 'site',
    geometry: 'ranges',
    target: 're/0',
    numbering: 'ordinal',
    value: { start: [2], stop: [4] },
  });
  const v = b.view({ ref: 'v', label: 'V', style: 'msa' });
  v.block({ alignment: 're' });
  v.track({ label: 'Rows' }).rows({});
  return b.getDocument();
}
const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const aln = (set: SeqSet, key: string) =>
  (set.byName(key, t) && t.kind === KeyKind.alignment ? set.alignments[t.id] : null) as Alignment;
const keys = (set: SeqSet) => set.mappings.map((m) => set.strings.get(m?.key ?? -1)).sort();

describe('computed alignments and mappings in the graph', () => {
  it("gives compile's set: the realigned rows, a mapping per row, the codon rows", async () => {
    const { ctx } = viewer();
    await ctx.setDocument(realignDoc());
    const graph = ctx.graph.assembled.value?.set as SeqSet;
    const reg = createRegistries();
    register(reg, [KernelRegistry, MvsRegistry]);
    const fresh = await compile(realignDoc(), reg, { load: () => Promise.reject(new Error('none')) });
    expect(aln(graph, 're').rowCount).toBe(3);
    expect(aln(graph, 're').columnCount).toBe(aln(fresh.set, 're').columnCount);
    expect(keys(graph)).toEqual(keys(fresh.set));
    expect(keys(graph)).toEqual(['cr.row-0', 'cr.row-1', 're.row-0', 're.row-1', 're.row-2']);
    expect(graph.annotations.some((a) => a && graph.strings.get(a.key) === 'site')).toBe(true);
    expect(ctx.graph.assembled.value?.view?.resolved.blocks.count).toBe(1);
  });
  it('recomputes the alignment when an input changes, and leaves it when something else does', async () => {
    const { ctx } = viewer();
    await ctx.setDocument(realignDoc());
    const before = aln(ctx.graph.assembled.value?.set as SeqSet, 're');
    await ctx.applyPatch(patch([{ op: 'set', ref: 'site', params: { value: { start: [3], stop: [5] } } }]));
    expect(aln(ctx.graph.assembled.value?.set as SeqSet, 're')).toBe(before);
    expect(refsRan(ctx)).not.toContain('re');
    await ctx.applyPatch(patch([{ op: 'set', ref: 's2', params: { value: 'VLAGWTQRSTEDNLK' } }]));
    expect(refsRan(ctx)).toContain('re');
    const after = aln(ctx.graph.assembled.value?.set as SeqSet, 're');
    expect(after).not.toBe(before);
    expect(after.columnCount).toBe(17);
    expect(after.rows.colToPos(2 as never, 2 as never)).toBe(0);
  });
});

describe('a query that names another query', () => {
  it("reads the set with the other's answer in it: no issue, however the nodes are ordered", async () => {
    // Regions of p53-like letters: the second holds three cysteines, the first one.
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: 'p', value: 'MAAAAAACAAACAAACAAACAAAA' });
    b.data.annotation({
      ref: 'regions',
      geometry: 'ranges',
      target: p,
      value: [
        { start: 0, stop: 8, label: 'one' },
        { start: 8, stop: 24, label: 'two' },
      ],
    });
    // The dependent one first in the document, so the order is the dependency's, not the document's.
    b.data.computed({
      ref: 'rich',
      provider: 'query',
      query: Q.read('(pick (items :of regions) :where (>= (count cysteines) 3))'),
    } as never);
    b.data.computed({ ref: 'cysteines', provider: 'query', query: Q.read('(and p (letters C))') } as never);
    b.view({ ref: 'v', label: 'V' }).track({ label: 'Letters' }).letters();
    const { ctx } = viewer();
    await ctx.setDocument(b.getDocument());
    const set = ctx.graph.assembled.value?.set as SeqSet;
    expect(ctx.issues.current.value.map((i) => i.message)).toEqual([]);
    const rich = set.annotations.find((a) => a && set.strings.get(a.key) === 'rich');
    expect([rich?.count, rich ? Array.from(rich.start) : null]).toEqual([1, [8]]);
  });
});

describe("edges and keys from core's `nodeRefs` (VS-24)", () => {
  const cysStarts = (ctx: SeqViewerContext) => {
    const s = ctx.graph.assembled.value?.set as SeqSet;
    const a = s.annotations.find((x) => x && s.strings.get(x.key) === 'cys');
    return a ? Array.from(a.start) : null;
  };
  const queryDoc = (p: string) => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: 'p', value: p });
    b.data.computed({ ref: 'cys', provider: 'query', query: Q.read('(and p (letters C))') });
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    return b.getDocument();
  };

  it('a query reruns when the sequence it reads changes, as a fresh open gives', async () => {
    const { ctx } = viewer();
    await ctx.setDocument(queryDoc('MACAAAAC'));
    expect(cysStarts(ctx)).toEqual([2, 7]);
    await ctx.applyPatch(patch([{ op: 'set', ref: 'p', params: { value: 'MCCCAAAA' } }]));
    const fresh = viewer();
    await fresh.ctx.setDocument(queryDoc('MCCCAAAA'));
    expect(cysStarts(ctx)).toEqual(cysStarts(fresh.ctx));
    expect(refsRan(ctx)).toContain('cys');
  });

  it('a `rules` theme over a computed column is made again when that column reruns', async () => {
    const b = createSeqViewSpec();
    b.data.alignment({ ref: 'msa', value: ['MKV-A', 'MKVLA', 'MRV-A'], alphabet: 'protein' } as never);
    b.data.computed({ ref: 'cons', provider: 'column_stats', alignment: 'msa', measure: 'conservation' } as never);
    const v = b.view({ ref: 'v', label: 'V', style: 'msa' } as never);
    v.block({ alignment: 'msa' } as never);
    const rules = [{ when: [{ measure: { column: 'cons' }, min: 0.6 }], color: '#ff0000' }];
    v.track({ ref: 'v.t' }).letters({ ref: 'v.l', color: { theme: 'rules', params: { rules } } } as never);
    const w = viewer();
    expect((await w.ctx.setDocument(b.getDocument())).ok).toBe(true);
    const made = w.created.filter((c) => c === 'rules:v.l').length;
    await w.ctx.applyPatch(patch([{ op: 'set', ref: 'cons', params: { measure: 'gap_fraction' } } as never]));
    expect(refsRan(w.ctx)).toContain('cons');
    expect(w.created.filter((c) => c === 'rules:v.l').length).toBe(made + 1);
  });
});
