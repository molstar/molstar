/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A headless export of a ruler in a numbering read from a file, on letters read from a file (_spec/README.md "Render",
// Export; Numberings from a source): the export opens the numbering's source although no layer's data names it, so the
// ruler draws its labels. Before, the source was never opened (the ruler's want waits for an open source) and the
// export said it settled with no tick drawn.

import { categoryColumn } from '@molstar/core/repr/tree';
import { encodeSeq, FastaRuntime, inlineSeq, SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { viewSources } from '../export.js';
import { renderSvg } from '../render-svg.js';

const Letters = 'EVQLVESGGGLVQPGGSLRL';
// biome-ignore format: a label a position, in index order
const Kabat = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '11', '12', '13', '14', '15', '52', '52A', '53', '54', '55', '56'];

const files: Record<string, string> = {
  'h.fa': `>h1 heavy\n${Letters}\n`,
  'kabat.seq': encodeSeq(
    inlineSeq(
      { name: 'h1', length: Letters.length, alphabet: 'protein' },
      {
        ranges: {
          form: 'list',
          start: Int32Array.from(Kabat, (_, i) => i),
          end: Int32Array.from(Kabat, (_, i) => i + 1),
          fields: [{ name: 'label', column: categoryColumn(Kabat) }],
        },
      },
    ),
  ),
};

const assets = {
  readText: async (uri: string) => files[uri],
  readBytes: async (uri: string) => new TextEncoder().encode(files[uri]),
};

function spec() {
  const b = createSeqSpec({ title: 'file-read numbering' });
  const fa = b.data.source('fa', { format: 'fasta', uri: 'h.fa' });
  const labels = b.data.source('labels', { format: 'seq', uri: 'kabat.seq' });
  const h = b.data.sequence('h', { alphabet: 'protein', source: fa, name: 'h1' }, (s) =>
    s.numbering('kabat', { method: 'labels', source: labels, scheme: 'kabat' }),
  );
  b.views.view('v', (v) => {
    const s = v.section('s', { of: h });
    s.track('ruler').layer('r', { type: 'ruler', numbering: 'kabat' });
    s.track('letters').layer('l', { type: 'letters' });
  });
  return b.build();
}

describe('a ruler in a file-read numbering, exported headless', () => {
  it('the view’s sources hold the numbering’s', () => {
    const model = resolveSpec(spec(), { runtimes: [FastaRuntime, SeqRuntime] });
    expect([...viewSources(model, resolveViews(model)[0])].sort()).toEqual(['/data/fa', '/data/labels']);
  });

  it('draws its labels, insertion codes included', async () => {
    const r = await renderSvg(spec(), {
      size: { width: 1400, height: 160 },
      runtimes: [FastaRuntime, SeqRuntime],
      assets,
      legend: 'none',
    });
    expect(r.settled).toBe(true);
    expect(r.issues.map((i) => i.message)).toEqual([]);
    const texts = [...r.svg.matchAll(/<text[^>]*>([^<]*)<\/text>/g)].map((m) => m[1]);
    for (const l of ['11', '52', '52A', '56']) expect(texts).toContain(l);
  });
});
