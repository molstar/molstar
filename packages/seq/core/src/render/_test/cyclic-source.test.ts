/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A cyclic sequence whose letters come from a source (_spec/README.md "Resolve", Systems): once the source is open its
// coverage's length is the model's, so a ring and the linear detail following it across the origin, inline rows that
// wrap the origin and a GFF3 row that ends past the length are drawn as they are over the same letters inline.

import type { Node } from '@molstar/core/repr/tree';
import { type FormatRuntime, inlineSeq, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { FakeClock } from '../../data/_test/synthetic.js';
import { resolveSpec, sourceLengths } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { type RangesLay, recordDrawLists, renderSvg } from '../index.js';
import { settle } from './helpers.js';

const N = 120;
const Letters = 'ACGTTGCA'.repeat(N / 8);
const Gff = [
  '##gff-version 3',
  'ring\tt\tCDS\t111\t130\t.\t+\t.\tID=across',
  'ring\tt\tCDS\t21\t30\t.\t+\t.\tID=inside',
].join('\n');

/** The spec over inline letters or a `seq` source of them (`fileLetters`: the file's, an edited source's). */
function spec(fromSource: boolean, fileLetters = Letters, format = 'seq'): Node {
  const b = createSeqSpec();
  const file = b.data.source('file', {
    format,
    value: seqToValue(
      inlineSeq(
        { name: 'ring', length: fileLetters.length, alphabet: 'dna', cyclic: true },
        { letters: { form: 'letters', start: 0, letters: fileLetters } },
      ),
    ),
  });
  const seq = b.data.sequence(
    'ring',
    fromSource
      ? { alphabet: 'dna', source: file, cyclic: true }
      : { alphabet: 'dna', letters: Letters, cyclic: true, name: 'ring' },
  );
  const sites = b.data.ranges('sites', {
    target: seq,
    rows: { start: [110, 40], end: [10, 50], id: ['wrap', 'plain'] },
  });
  const gff = b.data.source('gff', { format: 'gff3', value: Gff, fields: { category: 'type' } });
  const features = b.data.ranges('features', { target: seq, source: gff });
  b.views.view('main', (v) => {
    const ring = v.section('ring', { of: seq, layout: 'rings' });
    ring.track('sites').layer('sites', { type: 'ranges', data: sites });
    ring.track('features').layer('features', { type: 'ranges', data: features });
    ring.track('ruler').layer('ruler', { type: 'ruler' });
    const detail = v.section('detail', { of: seq, follow: ring });
    detail.track('sites').layer('sites', { type: 'ranges', data: sites });
    detail.track('features').layer('features', { type: 'ranges', data: features });
    detail.track('letters').layer('letters', { type: 'letters' });
  });
  return b.build();
}

const Runtimes = [SeqRuntime, RecordsRuntime];
const Size = { width: 720, height: 1000 };
// the ring's window over the origin: the detail follows it in two pieces
const Viewports = { '/views/main/ring': { start: 100, end: 140 } };

const pieces = (lay: RangesLay) =>
  Array.from({ length: lay.items }, (_, i) => `${lay.keys[i]}:${lay.starts[i]}-${lay.ends[i]}`).sort();

describe('a cyclic sequence read from a source', () => {
  it('resolves with the length its coverage gave: a wrapped row is two rows of one item, as inline', async () => {
    const before = resolveSpec(spec(true), { runtimes: Runtimes });
    expect(before.systems.get('/data/ring')?.length).toBeUndefined();
    expect(before.issues.map((i) => i.code)).toContain('seq.resolve.length.unknown');
    const s = await settle(spec(true), Runtimes);
    expect(sourceLengths(before, (src) => s.store.coverage(src)).size).toBe(1);
    expect(s.model.systems.get('/data/ring')?.length).toBe(N);
    expect(s.model.issues).toEqual([]);
    expect(s.model.annotations.get('/data/sites')?.items?.map((i) => `${i.key}:${i.start}-${i.end}`)).toEqual([
      'wrap:110-120',
      'wrap:0-10',
      'plain:40-50',
    ]);
    // the lengths are part of the model's hash
    expect(s.model.hash).not.toBe(before.hash);
  });

  it('a length is kept by its source: once the source is edited, the spec resolves with no length until it opens', async () => {
    const s = await settle(spec(true), Runtimes);
    const lengths = sourceLengths(s.model, (src) => s.store.coverage(src));
    expect(resolveSpec(spec(true), { runtimes: Runtimes, lengths }).systems.get('/data/ring')?.length).toBe(N);
    // another file under the same address: the old file's length is not taken
    const edited = spec(true, 'ACGT'.repeat(10));
    expect(resolveSpec(edited, { runtimes: Runtimes, lengths }).systems.get('/data/ring')?.length).toBeUndefined();
    expect((await settle(edited, Runtimes)).model.systems.get('/data/ring')?.length).toBe(40);
  });

  it('the ring, the detail across the origin and the GFF3 row past the length are drawn as over inline letters', async () => {
    const inline = await settle(spec(false), Runtimes);
    const sourced = await settle(spec(true), Runtimes);
    const a = await inline.frame({ size: Size, settled: true, viewports: Viewports });
    const b = await sourced.frame({ size: Size, settled: true, viewports: Viewports });
    expect(b.layout.issues).toEqual([]);
    const detail = b.layout.sections[1];
    expect(detail.axis?.pieces).toEqual([100, 120, 0, 20]);
    expect(pieces(detail.tracks[1].layers[0].lay as RangesLay)).toEqual(
      expect.arrayContaining(['across:0-10', 'across:110-120']),
    );
    expect(recordDrawLists(b.lists)).toBe(recordDrawLists(a.lists));
  });

  it('the export reads the source first and draws the same', async () => {
    const a = await renderSvg(spec(false), { size: Size, runtimes: Runtimes, viewports: Viewports });
    const b = await renderSvg(spec(true), { size: Size, runtimes: Runtimes, viewports: Viewports });
    expect(b.issues).toEqual([]);
    expect(b.settled).toBe(true);
    expect(recordDrawLists(b.lists)).toBe(recordDrawLists(a.lists));
  });

  it('a source that opens after the deadline never swaps the model under the render', async () => {
    // a `seq` source that opens only when released, after the export has given up on it
    let release: () => void = () => {};
    const late: FormatRuntime = {
      name: 'late',
      formats: ['late'],
      open: (entry, ctx) =>
        new Promise((resolve) => {
          release = () => resolve(SeqRuntime.open({ ...entry, format: 'seq' }, ctx));
        }),
    };
    const never: FormatRuntime = { name: 'never', formats: ['late'], open: () => new Promise(() => {}) };
    const render = async (runtime: FormatRuntime, after?: () => void) => {
      const clock = new FakeClock();
      const out = renderSvg(spec(true, Letters, 'late'), {
        size: Size,
        runtimes: [runtime, RecordsRuntime],
        viewports: Viewports,
        clock,
        timeoutMs: 100,
      });
      await clock.advance(100);
      after?.();
      await clock.advance(1000);
      return out;
    };
    const a = await render(late, () => release());
    const b = await render(never);
    expect(a.settled).toBe(false);
    expect(a.issues.map((i) => i.code)).toContain('seq.resolve.length.unknown');
    expect(a.svg).toBe(b.svg);
  });
});
