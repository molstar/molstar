/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Coverage } from '@molstar/seqx-core/data';
import { channelAddress } from '@molstar/seqx-core/data/channels';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { resolveSpec, sourceLengths } from '../model.js';
import { StatesSystemName, statesSystemOf } from '../model-region.js';

// Pairs, grids and matrices read by region (_spec/README.md "Resolve", Access plans by source): an access plan on the
// first system by its name, a field read as a channel; a matrix's states a system, its length from the inline states or
// its source's coverage (EC-220).

function spec() {
  const b = createSeqSpec();
  const src = b.data.source('hic', { format: 'seq', uri: 'hic.seq' });
  const vcf = b.data.source('vcf', { format: 'seq', uri: 'vcf.seq' });
  const chr1 = b.data.sequence('chr1', { alphabet: 'dna', source: src, name: 'chr1' });
  const chr2 = b.data.sequence('chr2', { alphabet: 'dna', source: src, name: 'chr2' });
  b.data.pairs('links', { a: chr1, b: chr2, source: src, field: 'score' } as never);
  b.data.grid('contacts', { x: chr1, source: src, field: 'balanced' } as never);
  b.data.matrix('genotypes', { target: chr1, source: vcf } as never);
  b.data.matrix('inline', {
    target: chr2,
    codes: new Uint8Array([0, 1, 2, 2, 1, 0]),
    states: ['A', 'B', 'C'],
  } as never);
  return b.build();
}

describe('access plans by source', () => {
  it('read pairs and a grid by region on their first system, a field as a channel (EC-208, EC-223)', () => {
    const m = resolveSpec(spec());
    const links = m.annotations.get('/data/links')!;
    expect(links.shape).toBe('pairs');
    expect(links.ySystem).toBe('/data/chr2');
    expect(links.pairs).toBeUndefined();
    expect(links.access).toEqual({
      kind: 'source',
      source: channelAddress('/data/hic', 'score'),
      system: 'chr1',
      shape: 'pairs',
    });
    const grid = m.annotations.get('/data/contacts')!;
    expect(grid.ySystem).toBe('/data/chr1');
    expect(grid.access).toMatchObject({ kind: 'source', system: 'chr1', shape: 'grid' });
  });

  it('give a matrix a states system, pending until its source’s coverage lists them (EC-220)', () => {
    const spec0 = spec();
    const m = resolveSpec(spec0);
    const g = m.annotations.get('/data/genotypes')!;
    expect(g.shape).toBe('matrix');
    expect(g.access).toMatchObject({ kind: 'source', source: '/data/vcf', system: 'chr1', shape: 'matrix' });
    const states = m.systems.get(statesSystemOf('/data/genotypes'))!;
    expect(states.name).toBe(StatesSystemName);
    expect(states.length).toBeUndefined();
    const coverage: Coverage = {
      systems: [{ name: 'chr1', length: 1000 }],
      shapes: [{ system: 'chr1', shape: 'matrix', states: ['NA12878', 'NA12891', 'NA12892'] }],
      levels: [1],
    };
    const lengths = sourceLengths(m, (s) => (s === '/data/vcf' ? coverage : undefined));
    const m2 = resolveSpec(spec0, { lengths });
    expect(m2.systems.get(statesSystemOf('/data/genotypes'))!.length).toBe(3);
  });

  it('hold an inline matrix as one item a position, its states inline', () => {
    const m = resolveSpec(spec());
    const a = m.annotations.get('/data/inline')!;
    expect(m.systems.get(statesSystemOf('/data/inline'))!.length).toBe(3);
    expect(a.access?.kind).toBe('inline');
    expect(a.access?.shape).toBe('matrix');
  });
});
