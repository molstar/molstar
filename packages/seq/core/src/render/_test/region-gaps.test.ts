/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn } from '@molstar/core/repr/tree';
import { type PairsList, type SeqData, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import type { MatrixLay } from '../matrix.js';
import { settle } from './helpers.js';

// Region reads' gaps the stories ran into (the lane of the region reads' stories found them; core's behaviour is not
// theirs to change): each a test of what a reader expects, marked `fails` while core does not do it, so the quick lane
// stays green and a fix turns the mark into a failure to remove.

beforeAll(async () => {
  await loadGenome();
});

describe('region reads: gaps', () => {
  it('a matrix of bins is drawn when its system ends inside the last tile (1,000,000 positions, tiles of 524,288)', async () => {
    // 1,000,000 positions at bins of 2,048 (524,288 a tile): the table is loaded and no image is made; at 524,288 or 1,048,576
    // positions, whole tiles, it is
    const states = 2;
    const bins = 489;
    const data: SeqData = {
      coverage: {
        systems: [{ name: 'chr', length: 1_000_000 }],
        shapes: [{ system: 'chr', shape: 'matrix', states: ['a', 'b'] }],
        levels: [1, 2048],
      },
      pairs: [
        {
          request: { system: 'chr', start: 0, end: 1_000_000, level: 1, shapes: ['matrix'] },
          answer: {
            ok: true,
            tables: {
              matrix: {
                form: 'items',
                start: Int32Array.from([10]),
                end: Int32Array.from([11]),
                states,
                codes: Uint8Array.from([1, 2]),
              },
            },
          },
        },
        {
          request: { system: 'chr', start: 0, end: 1_000_000, level: 2048, shapes: ['matrix'] },
          answer: {
            ok: true,
            tables: {
              matrix: {
                form: 'summary',
                start: 0,
                states,
                classes: new Uint8Array(bins * states).fill(1),
                shares: new Uint8Array(bins * states).fill(255),
              },
            },
          },
        },
      ],
    };
    const b = createSeqSpec({ title: 'partial bin' });
    const src = b.data.source('src', { format: 'seq', value: seqToValue(data) as never });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
    const m = b.data.matrix('m', { target: chr, source: src } as never);
    b.views.view('v', (v) => {
      v.section('chr', { of: chr }).track('m', { height: 60 }).layer('m', { type: 'matrix', data: m });
    });
    const s = await settle(b.build(), [SeqRuntime]);
    const f = await s.frame({ size: { width: 720, height: 200 } });
    const lay = f.layout.sections[0].tracks[0].layers[0].lay as unknown as MatrixLay;
    expect(lay.level).toBe(2048);
    expect(lay.images.length).toBe(1);
  });

  it('a BEDPE source is read by a pairs layer: its pairs shape listed, its bins by its own levels', async () => {
    const lines = Array.from({ length: 20 }, (_, i) =>
      ['chr', 1000 * i, 1000 * i + 10, 'chr', 90_000 - 1000 * i, 90_010 - 1000 * i, `l${i}`, 1, '+', '+'].join('\t'),
    );
    const b = createSeqSpec({ title: 'bedpe' });
    const sizes = b.data.source('sizes', { format: 'chrom-sizes', value: 'chr\t2000000\n' as never });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: sizes, name: 'chr' });
    const bp = b.data.source('bp', { format: 'bedpe', value: `${lines.join('\n')}\n` as never });
    const pairs = b.data.pairs('p', { a: chr, source: bp } as never);
    b.views.view('v', (v) => {
      v.section('chr', { of: chr }).track('t').layer('p', { type: 'pairs', data: pairs });
    });
    const s = await settle(b.build(), [BedRuntime]);
    const f = await s.frame({ size: { width: 720, height: 300 } });
    const lay = f.layout.sections[0].sectionLays[0].lay as unknown as { n: number };
    // the `bed` runtime lists `pairs` and levels of its own (its bin pairs per tile), so the view at the fit reads a
    // coarse level within a want's tiles, never 1,954 tiles at level 1 (`seq.store.too_many_tiles`)
    expect(s.store.issues.value.map((i) => i.code)).not.toContain('seq.store.too_many_tiles');
    expect(lay.n).toBe(20);
  });

  it('pairs read by region over an assembly are laid on its blocks (a ring’s chords between chromosomes)', async () => {
    const list: PairsList = {
      form: 'list',
      aSystem: categoryColumn(['c1', 'c1']),
      bSystem: categoryColumn(['c2', 'c2']),
      aStart: Int32Array.from([10, 40]),
      aEnd: Int32Array.from([20, 50]),
      bStart: Int32Array.from([30, 60]),
      bEnd: Int32Array.from([40, 70]),
    };
    const data: SeqData = {
      coverage: {
        systems: [
          { name: 'c1', length: 100 },
          { name: 'c2', length: 100 },
        ],
        shapes: [
          { system: 'c1', shape: 'pairs' },
          { system: 'c2', shape: 'pairs' },
        ],
        levels: [1],
      },
      pairs: ['c1', 'c2'].map((system) => ({
        request: { system, start: 0, end: 100, level: 1, shapes: ['pairs' as const] },
        answer: { ok: true as const, tables: { pairs: list } },
      })),
    };
    const b = createSeqSpec({ title: 'pairs over an assembly' });
    const g = b.data.assembly(
      'toy',
      assembly('Toy', [
        { name: 'c1', length: 100 },
        { name: 'c2', length: 100 },
      ]),
    );
    const src = b.data.source('src', { format: 'seq', value: seqToValue(data) as never });
    const pairs = b.data.pairs('links', { a: g, source: src } as never);
    b.views.view('v', (v) => {
      genome(v, 'genome', g, { layout: 'rings' }, (s) => {
        s.track('links').layer('links', { type: 'pairs', data: pairs });
      });
    });
    const s = await settle(b.build(), [SeqRuntime]);
    const f = await s.frame({ size: { width: 720, height: 600 } });
    expect(f.layout.sections[0].sectionLays.length).toBe(1);
    // the two links, each a chord from c1 to c2 (found from either end, once)
    const lay = f.layout.sections[0].sectionLays[0].lay as unknown as { n: number };
    expect(lay.n).toBe(2);
  });

  it('held pairs past the gate bin at the pairs level a source is read at: held and from a file one picture (RS gap 6)', async () => {
    // 6,000 pairs between 60 hotspots of a 2 Mb chromosome, held in the spec and served at level 1 by a `seq` source
    const n = 6000;
    let x = 7;
    const r = () => (x = (x * 1103515245 + 12345) % 2 ** 31) / 2 ** 31;
    const as = Int32Array.from({ length: n }, () => Math.floor(r() * 60) * 33_000 + Math.floor(r() * 500));
    const bs = Int32Array.from({ length: n }, () => Math.floor(r() * 60) * 33_000 + Math.floor(r() * 500));
    const list: PairsList = {
      form: 'list',
      aStart: as,
      aEnd: as.map((v) => v + 10),
      bStart: bs,
      bEnd: bs.map((v) => v + 10),
    };
    const data: SeqData = {
      coverage: {
        systems: [{ name: 'chr', length: 2_000_000 }],
        shapes: [{ system: 'chr', shape: 'pairs' }],
        levels: [1],
      },
      pairs: [
        {
          request: { system: 'chr', start: 0, end: 2_000_000, level: 1, shapes: ['pairs'] },
          answer: { ok: true, tables: { pairs: list } },
        },
      ],
    };
    const b = createSeqSpec({ title: 'held and read' });
    const sizes = b.data.source('sizes', { format: 'chrom-sizes', value: 'chr\t2000000\n' as never });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: sizes, name: 'chr' });
    const held = b.data.pairs('held', {
      a: chr,
      rows: { a_start: list.aStart, a_end: list.aEnd, b_start: list.bStart, b_end: list.bEnd },
    } as never);
    const src = b.data.source('src', { format: 'seq', value: seqToValue(data) as never });
    const read = b.data.pairs('read', { a: chr, source: src } as never);
    b.views.view('v', (v) => {
      const s = v.section('chr', { of: chr });
      s.track('held').layer('held', { type: 'pairs', data: held });
      s.track('read').layer('read', { type: 'pairs', data: read });
    });
    const s = await settle(b.build(), [SeqRuntime, BedRuntime]);
    let f = await s.frame({ size: { width: 720, height: 300 } });
    await s.store.whenSettled();
    f = await s.frame({ size: { width: 720, height: 300 } });
    type Lay = {
      prep?: { level?: number; weight?: Float32Array; data?: { n: number; aStart: Int32Array; bStart: Int32Array } };
    };
    const [h, d] = f.layout.sections[0].sectionLays.map((l) => l.lay as unknown as Lay);
    const bins = (l: Lay) => {
      const out = new Map<string, number>();
      for (let k = 0; k < (l.prep?.data?.n ?? 0); k++) {
        const key = `${l.prep!.data!.aStart[k]}|${l.prep!.data!.bStart[k]}`;
        out.set(key, (out.get(key) ?? 0) + (l.prep!.weight?.[k] ?? 1));
      }
      return out;
    };
    expect(d.prep?.weight).toBeDefined();
    expect(h.prep?.weight).toBeDefined();
    expect(h.prep?.level).toBe(d.prep?.level);
    expect(bins(h)).toEqual(bins(d));
  });
});
