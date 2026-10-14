/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST7: streamed data placed through a path (seqspec:sources#data-patches): a genome's genes reach a protein's block by
// its codon map, but not the copies of a merged block each reached by a path of its own (`svs.layer.streamed_merged`):
// their regions would land on one set of indices.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { compile, resolveView, supportedAll } from '../../src/format/index.js';
import { registries } from './load-fixtures.js';
import { CODING_CDS, codingProtein, GRCH38, hg38Streams } from './stream-provider.js';

const reg = registries();

/** P1 twice (units A and B), each coded by its own CDS on chr2, its genes shown on one section of `copies`. */
function copiesDoc(copies: 'one' | 'each', sameCds = false) {
  const b = createSeqViewSpec({ title: 'Copies' });
  const genome = b.data.source({ ref: 'genome', stream: 'genome', data_kind: 'letters' });
  b.data.complex({ ref: 'hg38', source_ref: genome, label: 'GRCh38' });
  b.data.source({ ref: 'genes_src', stream: 'genes', data_kind: 'items', coordinate_system: 'hg38' });
  b.data.annotation({ ref: 'genes', geometry: 'ranges', source_ref: 'genes_src' });
  const p = b.data.sequence({ ref: 'p1', alphabet: 'protein', value: codingProtein() });
  const a = b.data.unit({ ref: 'p1.A', sequence: p, label: 'A' });
  const u = b.data.unit({ ref: 'p1.B', sequence: p, label: 'B' });
  const shifted = CODING_CDS.map(([s, t]): [number, number] => [s + 30, t + 30]);
  b.data.mapping({ ref: 'a.codon', from: a, to: 'hg38/chr2', method: 'codon', cds: CODING_CDS });
  b.data.mapping({ ref: 'b.codon', from: u, to: 'hg38/chr2', method: 'codon', cds: sameCds ? CODING_CDS : shifted });
  const s = b.view({ ref: 'v', label: 'Copies' }).section({ ref: 'v.main', label: 'P1', copies });
  s.block({ unit: a });
  s.block({ unit: u });
  s.track({ label: 'Genes' }).ranges({ data: 'genes' });
  return b.getDocument();
}

async function resolved(doc: ReturnType<typeof copiesDoc>) {
  const s = hg38Streams({ chromosomes: GRCH38.slice(0, 3) });
  const loaders = { load: () => Promise.reject(new Error('none')), stream: s.stream };
  const { set, issues } = await compile(doc, reg, loaders);
  expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
  const r = resolveView(doc, set, 'v', reg, supportedAll(reg));
  return { codes: r.issues.map((i) => i.code), placements: r.view.layers.placementCount[0] };
}

describe('streamed data through paths onto copies', () => {
  it('each copy shown apart: each block through its own path', async () => {
    expect(await resolved(copiesDoc('each'))).toEqual({ codes: [], placements: 2 });
  });

  it('the copies merged, each by its own path: not placed, and said', async () => {
    expect(await resolved(copiesDoc('one'))).toEqual({ codes: ['svs.layer.streamed_merged'], placements: 0 });
  });
});
