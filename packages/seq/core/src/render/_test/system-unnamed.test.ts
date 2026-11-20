/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A sequence read from a source whose system is not named (_spec/README.md "Resolve", Systems: A length from coverage):
// once the source is open and names no system after the sequence, its length stays unknown and nothing over it lays
// out; `seq.resolve.system.unnamed` says so, an error (a spec fault, never a wait), and a sequence named after the
// source's system has none.

import type { Node } from '@molstar/core/repr/tree';
import { DataStore, inlineSeq, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec, sourceLengths } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { openModel } from '../export.js';

const Letters = 'ACGT'.repeat(30);

/** A section over a sequence read from a `seq` source whose one system is `chr1`; the sequence named `name`. */
function spec(name: string | undefined, id = 'chr1'): Node {
  const b = createSeqSpec();
  const file = b.data.source('file', {
    format: 'seq',
    value: seqToValue(
      inlineSeq(
        { name: 'chr1', length: Letters.length, alphabet: 'dna' },
        { letters: { form: 'letters', start: 0, letters: Letters } },
      ),
    ),
  });
  const seq = b.data.sequence(id, { alphabet: 'dna', source: file, ...(name ? { name } : {}) });
  b.views.view('v', (v) => v.section('s', { of: seq }).track('l').layer('l', { type: 'letters' }));
  return b.build();
}

const codes = (issues: readonly { code: string }[]) => issues.map((i) => i.code);

describe('a sequence whose source names no system after it', () => {
  it('is said once the source is open (an error at the sequence), never before', async () => {
    const before = resolveSpec(spec(undefined, 'genome'), { runtimes: [SeqRuntime] });
    expect(codes(before.issues)).not.toContain('seq.resolve.system.unnamed');
    const opened = await openModel(spec(undefined, 'genome'), new DataStore(), { runtimes: [SeqRuntime] });
    const unnamed = opened.issues.filter((i) => i.code === 'seq.resolve.system.unnamed');
    expect(unnamed).toHaveLength(1);
    expect(unnamed[0]).toMatchObject({ severity: 'error', address: '/data/genome', path: ['params', 'source'] });
    expect(unnamed[0].message).toContain('genome');
    expect(opened.model.systems.get('/data/genome')?.length).toBeUndefined();
  });

  it('a `name` that is not the source’s system is the param said; one that is has no issue and a length', async () => {
    const wrong = await openModel(spec('chr2', 'genome'), new DataStore(), { runtimes: [SeqRuntime] });
    expect(wrong.issues.find((i) => i.code === 'seq.resolve.system.unnamed')?.path).toEqual(['params', 'name']);
    const right = await openModel(spec('chr1', 'genome'), new DataStore(), { runtimes: [SeqRuntime] });
    expect(codes(right.issues)).toEqual([]);
    expect(right.model.systems.get('/data/genome')?.length).toBe(Letters.length);
    const byId = await openModel(spec(undefined), new DataStore(), { runtimes: [SeqRuntime] });
    expect(codes(byId.issues)).toEqual([]);
  });

  it('a source not open yet, or naming no systems, says nothing', () => {
    const model = resolveSpec(spec(undefined, 'genome'), { runtimes: [SeqRuntime] });
    expect(sourceLengths(model, () => undefined).size).toBe(0);
    expect(sourceLengths(model, () => ({ systems: [], shapes: [], levels: [1] })).size).toBe(0);
  });
});
