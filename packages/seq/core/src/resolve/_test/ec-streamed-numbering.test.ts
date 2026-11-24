/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// EC-195: the default numbering of a streamed sequence (1-based labels) is stated
// by the spec and the model, known before any of its data arrives, not learned from a resolver.

import { checkSeqSpec, createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveAddressIn } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { resolveSpec } from '../model.js';

function streamed() {
  const b = createSeqSpec({ title: 'streamed' });
  const src = b.data.source('chr', { format: 'fasta', uri: 'x.fa' });
  b.data.sequence('chr1', { alphabet: 'dna', source: src });
  return b.build();
}

describe('EC-195: the default numbering of a streamed sequence', () => {
  const root = streamed();
  // no runtime and no data: nothing was read
  const model = resolveSpec(root);

  it('is stated by the model without data: the length open, `ordinal` the label numbering, `index` the 0-based', () => {
    expect(checkSeqSpec(root).filter((i) => i.severity === 'error')).toEqual([]);
    const sys = model.systems.get('/data/chr1')!;
    expect(sys.length).toBeUndefined();
    expect(sys.numberings).toEqual(['index', 'ordinal']);
    // the only thing said is that no runtime reads the source: nothing about its numbering
    expect(model.issues.map((i) => i.code)).toEqual(['seq.resolve.source.runtime']);
  });

  it('a label is read 1-based by default, 0-based only by the `index` numbering, before any data', () => {
    const ranges = (address: string) => {
      const r = resolveAddressIn(root, model, address);
      const loc = r.result?.type === 'location' ? r.result.location : undefined;
      return { ranges: loc ? Array.from(loc.parts[0]?.ranges ?? []) : undefined, codes: r.issues.map((i) => i.code) };
    };
    expect(ranges('/data/chr1/(at :label 1)')).toEqual({ ranges: [0, 1], codes: [] });
    expect(ranges('/data/chr1/(at :label 500)')).toEqual({ ranges: [499, 500], codes: [] });
    expect(ranges('/data/chr1/(at :label 500 :numbering ordinal)')).toEqual({ ranges: [499, 500], codes: [] });
    expect(ranges('/data/chr1/(at :label 500 :numbering index)')).toEqual({ ranges: [500, 501], codes: [] });
    expect(ranges('/data/chr1/(range :from 10 :to 20)').ranges).toEqual([9, 20]);
  });
});
