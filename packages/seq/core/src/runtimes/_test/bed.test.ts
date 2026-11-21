/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { gzipSync } from 'node:zlib';
import { fileURLToPath } from 'node:url';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { columnValues } from '@molstar/core/repr/tree';
import { bgzip } from '@molstar/test-harness/genome/bgzf-writer';
import { describe, expect, it, vi } from 'vitest';
import { sourceSuite } from '../../data/_test/source-suite.js';
import { delayedRuntime } from '../../data/async-source.js';
import { DefaultOpenLimits, type OpenContext, type OpenedSource, type SourceEntry } from '../../data/protocol.js';
import { SourceOpenError } from '../../data/seq-runtime.js';
import { checkAnswer } from '../../data/tables.js';
import { type BedpeSource, BedRuntime, bedRuntime } from '../bed.js';
import { indexPairs, PairsBuilder, pairsIn, pairsSchema } from '../bed-pairs.js';
import { BedCauses } from '../bed-issues.js';
import {
  generateBed6,
  generateBed12,
  generateBedGraph,
  generateBedpe,
  indexed,
  memoryAssets,
  plainRows,
} from './bed-fixtures.js';

// The BED family's runtime (_spec/README.md "Data", BED family): the source suite on generated files and on the
// recorded real slice (UCSC hg38 cytobands), whole and by region; a region read equals the whole read clipped (a
// property over generated files and random regions); fuzzed and truncated files without a throw; BED12's pieces,
// `fields` by number, identities; BEDPE found from either end, each pair once, its bytes bounded; chrom.sizes.

const signal = new AbortController().signal;
const names = ['chr1', 'chr2', 'chrX'];

function context(files: Record<string, Uint8Array | string>, ranges = true, limits = DefaultOpenLimits): OpenContext {
  return { assets: memoryAssets(files, ranges), signal, limits };
}

const entry = (format: string, more: Partial<SourceEntry> = {}): SourceEntry => ({
  address: '/data/s',
  format,
  uri: 'f',
  ...more,
});

/** The rejection's issue codes of an open, or 'opened'. */
async function opened(p: Promise<OpenedSource>): Promise<string[]> {
  try {
    const s = await p;
    s.dispose();
    return ['opened'];
  } catch (e) {
    if (e instanceof SourceOpenError) return e.issues.map((i) => i.code);
    throw e;
  }
}

/** The real slice: UCSC hg38's `cytoBandIdeo.txt.gz`, recorded as pages (examples' `recorded-pages`), read whole. */
const Recording = fileURLToPath(
  new URL('../../../../examples/fixtures/genomes/cytoBandIdeo.txt.gz.pages.json', import.meta.url),
);
function recordedBytes(): Uint8Array {
  const r = JSON.parse(readFileSync(Recording, 'utf8')) as {
    size: number;
    pageSize: number;
    pages: Record<string, string>;
  };
  const out = new Uint8Array(r.size);
  for (const [k, b64] of Object.entries(r.pages)) out.set(Buffer.from(b64, 'base64'), Number(k) * r.pageSize);
  return out;
}

describe('the source suite', () => {
  const bed = generateBed12(600, names, 4);
  const graph = generateBedGraph(500, names, 5);

  it('BED read whole: plain, gzip and BGZF by their magic, with and without range reads', async () => {
    const text = new TextEncoder().encode(bed.text);
    for (const bytes of [text, new Uint8Array(gzipSync(text)), bgzip(text, 3000).bytes])
      for (const ranges of [true, false]) await sourceSuite(BedRuntime, entry('bed'), context({ f: bytes }, ranges));
  });

  it('BED and bedGraph by region through tabix (level 1)', async () => {
    const b = indexed(bed);
    await sourceSuite(BedRuntime, entry('bed', { index: 'i' }), context({ f: b.file, i: b.index }));
    await sourceSuite(BedRuntime, entry('bed', { index: 'i' }), context({ f: b.file, i: b.index }, false));
    const g = indexed(graph);
    await sourceSuite(BedRuntime, entry('bedgraph', { index: 'i' }), context({ f: g.file, i: g.index }));
  });

  it('bedGraph read whole, and the synthetic twin delayed', async () => {
    await sourceSuite(BedRuntime, entry('bedgraph'), context({ f: graph.text }));
    await sourceSuite(delayedRuntime(BedRuntime, { delayMs: 1 }), entry('bed'), context({ f: bed.text }));
  });

  it('the recorded real slice: hg38’s cytobands, a band a row, `stain` by its column number (EC-190)', async () => {
    const bytes = recordedBytes();
    const source = entry('bed', { fields: { stain: '4' } });
    await sourceSuite(BedRuntime, source, context({ f: bytes }));
    const s = await BedRuntime.open(source, context({ f: bytes }));
    expect(s.issues).toEqual([]);
    expect(s.coverage.systems.length).toBeGreaterThan(24);
    expect(s.coverage.shapes[0].fields?.map((f) => f.name)).toEqual(['stain', 'label']);
    // one tile holds a whole chromosome: the level a navigator's ideogram reads
    const level = s.coverage.levels.at(-1)!;
    const a = await s.read({ system: 'chr17', start: 0, end: 1024 * level, level, shapes: ['ranges'] }, signal);
    expect(a.ok).toBe(true);
    if (!a.ok || a.tables.ranges?.form !== 'list') return;
    const r = a.tables.ranges;
    const stain = columnValues(r.fields!.find((f) => f.name === 'stain')!.column);
    expect(stain).toContain('acen');
    expect(stain).toContain('gneg');
    expect(columnValues(r.id!)[0]).toBe('p13.3@0-3400000');
    expect(r.end[r.end.length - 1]).toBe(83257441);
  });
});

describe('a region read equals the whole read clipped', () => {
  // positions to 400 kb (BED) and 200 kb (bedGraph): the regions are drawn inside the data
  for (const [format, kind, gen, span] of [
    ['bed', 'BED12', (seed: number) => generateBed12(400 + seed * 37, names, seed), 420_000],
    ['bed', 'BED6 with items of no length', (seed: number) => generateBed6(400 + seed * 37, names, seed), 420_000],
    ['bedgraph', 'bedGraph', (seed: number) => generateBedGraph(300 + seed * 41, names, seed), 200_000],
  ] as const) {
    it(`${kind}: generated files, random regions, every answer equal row for row`, async () => {
      const rng = new PCG(11);
      for (let seed = 1; seed <= 6; seed++) {
        const g = gen(seed);
        const { file, index } = indexed(g, 1000 + seed * 500);
        const whole = await BedRuntime.open(entry(format), context({ f: g.text }));
        const region = await BedRuntime.open(entry(format, { index: 'i' }), context({ f: file, i: index }));
        if (format === 'bed') expect(region.coverage.shapes[0].fields).toEqual(whole.coverage.shapes[0].fields);
        for (let k = 0; k < 25; k++) {
          const pick = names[Math.floor(rng.float() * names.length)];
          // every fifth region starts at a line's own start (an item of no length there included)
          const at = g.lines[Math.floor(rng.float() * g.lines.length)];
          const start = k % 5 === 1 ? at.begin : Math.floor(rng.float() * span);
          const end = start + Math.floor(rng.float() * (k % 5 === 0 ? 2 : 30_000));
          const system = k % 5 === 1 ? names[at.seq] : pick;
          const request = { system, start, end, level: 1, shapes: [format === 'bed' ? 'ranges' : 'values'] as const };
          const w = await whole.read({ ...request, shapes: [...request.shapes] }, signal);
          const r = await region.read({ ...request, shapes: [...request.shapes] }, signal);
          expect(checkAnswer(r, { ...request, shapes: [...request.shapes] })).toBeUndefined();
          expect(plainRows(r), `${format} ${seed} ${system}:${start}-${end}`).toEqual(plainRows(w));
        }
        whole.dispose();
        region.dispose();
      }
    });
  }

  it('the first data line read past `track`, `browser` and `#` lines where the index’s first chunk starts on them', async () => {
    const g = generateBed12(120, names, 31);
    const head = 'track name=t\nbrowser position chr1\n# c\n';
    const n = new TextEncoder().encode(head).length;
    // the header lines indexed as the data's first chunk's start (an index whose first line is a `track` line)
    const lines = [
      { seq: 0, begin: 0, end: 1, start: 0, stop: n },
      ...g.lines.map((l) => ({ ...l, start: l.start + n, stop: l.stop + n })),
    ];
    const { file, index } = indexed({ text: head + g.text, lines, names: g.names });
    const region = await BedRuntime.open(entry('bed', { index: 'i' }), context({ f: file, i: index }));
    const whole = await BedRuntime.open(entry('bed'), context({ f: g.text }));
    expect(region.coverage.shapes[0].fields).toEqual(whole.coverage.shapes[0].fields);
    expect(region.coverage.shapes[0].fields?.length).toBe(6);
  });
});

describe('fuzzed and truncated files read without a throw (EC-40, EC-92)', () => {
  const bed = generateBed12(200, names, 9);
  const text = new TextEncoder().encode(bed.text);

  it('bytes flipped, lines cut, files truncated: opened or refused with an issue, every answer checked', async () => {
    const rng = new PCG(21);
    const variants: Uint8Array[] = [];
    for (let k = 0; k < 40; k++) {
      const v = text.slice(0, Math.floor(rng.float() * text.length));
      for (let f = 0; f < 20; f++) v[Math.floor(rng.float() * v.length)] = Math.floor(rng.float() * 256);
      variants.push(v);
    }
    const gz = bgzip(text, 2000).bytes;
    for (let k = 0; k < 20; k++) {
      const v = gz.slice(0, Math.floor(rng.float() * gz.length));
      if (k % 2) v[Math.floor(rng.float() * v.length)] ^= 0x55;
      variants.push(v, new Uint8Array(gzipSync(text)).slice(0, Math.floor(rng.float() * 2000)));
    }
    for (const format of ['bed', 'bedgraph', 'bedpe', 'chrom-sizes'])
      for (const v of variants) {
        const got = await opened(BedRuntime.open(entry(format), context({ f: v })));
        expect(got.length).toBeGreaterThan(0);
        if (got[0] !== 'opened') continue;
        const s = await BedRuntime.open(entry(format), context({ f: v }));
        for (const sh of s.coverage.shapes.slice(0, 3)) {
          const request = { system: sh.system, start: 0, end: 50_000, level: 1, shapes: [sh.shape] };
          expect(checkAnswer(await s.read(request, signal), request)).toBeUndefined();
        }
      }
  });

  it('a truncated or damaged BGZF by region: the open refused or regions answered as failures, never wrong rows', async () => {
    const { file, index } = indexed(bed, 1500);
    const rng = new PCG(5);
    for (let k = 0; k < 15; k++) {
      const v = file.slice(0, file.length - 30 - Math.floor(rng.float() * file.length * 0.5));
      if (k % 3 === 0) v[Math.floor(rng.float() * v.length)] ^= 0xff;
      const got = await opened(BedRuntime.open(entry('bed', { index: 'i' }), context({ f: v, i: index })));
      if (got[0] !== 'opened') {
        expect(got[0]).toMatch(/^seq\.source\./);
        continue;
      }
      const s = await BedRuntime.open(entry('bed', { index: 'i' }), context({ f: v, i: index }));
      expect(s.issues?.map((i) => i.code)).toContain('seq.source.bgzf.eof');
      for (const system of names) {
        const request = { system, start: 0, end: 500_000, level: 1, shapes: ['ranges' as const] };
        const a = await s.read(request, signal);
        expect(checkAnswer(a, request)).toBeUndefined();
      }
    }
  });
});

describe('BED', () => {
  const read = async (text: string, more: Partial<SourceEntry> = {}, system = 'chr1', start = 0, end = 10_000) => {
    const s = await BedRuntime.open(entry('bed', more), context({ f: text }));
    const a = await s.read({ system, start, end, level: 1, shapes: ['ranges'] }, signal);
    return { s, a, t: a.ok && a.tables.ranges?.form === 'list' ? a.tables.ranges : undefined };
  };

  it('BED12’s blocks are pieces of one item, each with its extent; identity name@start-end (EC-196)', async () => {
    const line = 'chr1\t100\t1000\ttx1\t500\t-\t200\t900\t255,0,0\t3\t100,200,100,\t0,400,800,';
    const { t, s } = await read(`${line}\n`);
    // BED's name is the label (the items' label field), its strand a field too (a layer's filter reads it)
    expect(s.coverage.shapes[0].fields?.map((f) => `${f.name}:${f.type}`)).toEqual([
      'label:category',
      'score:float32',
      'thick_start:int32',
      'thick_end:int32',
      'item_rgb:category',
      'strand:int8',
    ]);
    expect(Array.from(t!.start)).toEqual([100, 500, 900]);
    expect(Array.from(t!.end)).toEqual([200, 700, 1000]);
    expect(Array.from(t!.itemStart!)).toEqual([100, 100, 100]);
    expect(Array.from(t!.itemEnd!)).toEqual([1000, 1000, 1000]);
    expect(Array.from(t!.strand!)).toEqual([-1, -1, -1]);
    expect(columnValues(t!.id!)).toEqual(['tx1@100-1000', 'tx1@100-1000', 'tx1@100-1000']);
    // a region in the intron between the first two pieces meets none of them
    const s2 = await BedRuntime.open(entry('bed'), context({ f: `${line}\n` }));
    const gap = await s2.read({ system: 'chr1', start: 250, end: 450, level: 1, shapes: ['ranges'] }, signal);
    expect(plainRows(gap)).toEqual([]);
  });

  it('blocks that do not fit their item leave the line out, said with its line number', async () => {
    const bad = [
      'chr1\t0\t100\ta\t0\t+\t0\t100\t0\t2\t10,10,\t0,',
      'chr1\t0\t100\tb\t0\t+\t0\t100\t0\t2\t10,10,\t50,0,',
      'chr1\t0\t100\tc\t0\t+\t0\t100\t0\t1\t200,\t0,',
      'chr1\t0\t100\td\t0\t+\t0\t100\t0\t1\t100,\t0,',
      // UCSC's rule: a block of no length, one not from the item's start, one short of its end, two overlapping
      'chr1\t0\t100\te\t0\t+\t0\t100\t0\t2\t0,100,\t0,0,',
      'chr1\t0\t100\tf\t0\t+\t0\t100\t0\t1\t90,\t10,',
      'chr1\t0\t100\tg\t0\t+\t0\t100\t0\t2\t10,10,\t0,50,',
      'chr1\t0\t100\th\t0\t+\t0\t100\t0\t2\t60,50,\t0,50,',
    ];
    const { t, s } = await read(`${bad.join('\n')}\n`);
    expect(columnValues(t!.id!)).toEqual(['d@0-100']);
    const lines = s.issues?.filter((i) => i.code === 'seq.source.bed.line');
    expect(lines?.map((i) => i.message)).toEqual([
      'Blocks whose counts are not `block_count`: 1 (line 1)',
      'Blocks out of order, overlapping or outside the item: 4 (line 2, line 3, line 6, line 8)',
      'A block of no length: 1 (line 5)',
      'Blocks that do not reach the item’s end: 1 (line 7)',
    ]);
    expect(lines?.[0].severity).toBe('warning');
  });

  it('`bed N+M` by `fields`, a column by its 0-based number; the column not read under BED’s name (EC-190)', async () => {
    const text = 'chr1\t0\t2300000\tp36.33\tgneg\nchr1\t2300000\t5300000\tp36.32\tgpos25\n';
    const plain = await read(text);
    // without the mapping, column 4 is BED's score: `gneg` is no number, said
    expect(plain.s.issues?.map((i) => i.code)).toEqual(['seq.source.bed.cells']);
    const { t, s } = await read(text, { fields: { stain: '4' } }, 'chr1', 0, 6_000_000);
    expect(s.issues).toEqual([]);
    expect(t!.fields!.map((f) => f.name)).toEqual(['stain', 'label']);
    expect(columnValues(t!.fields![0].column)).toEqual(['gneg', 'gpos25']);
    // narrowPeak (`bed 6+4`): signal and p by number, typed by `types`; BED's names stop at the first mapped column, so
    // columns 8 and 9 are `col8` and `col9`, never `item_rgb` and `block_count`; a column past the first line's refused
    const peak = 'chr1\t10\t20\tp1\t0\t.\t5.5\t3.2\t-1\t4\n';
    const n = await read(peak, { fields: { signal: '6', p: '7' }, types: { signal: 'float32', p: 'float32' } });
    expect(n.t!.fields!.map((f) => `${f.name}:${columnValues(f.column)[0]}`)).toEqual([
      'signal:5.5',
      'p:3.200000047683716',
      'label:p1',
      'score:0',
      'col8:-1',
      'col9:4',
      'strand:0',
    ]);
    expect(await opened(BedRuntime.open(entry('bed', { fields: { x: '12' } }), context({ f: peak })))).toEqual([
      'seq.source.bed.field',
    ]);
    expect(await opened(BedRuntime.open(entry('bed', { fields: { x: 'nope' } }), context({ f: peak })))).toEqual([
      'seq.source.bed.field',
    ]);
  });

  it('a format’s base is fixed: `base` and `ends` not read, said once, whole and by region (EC-195)', async () => {
    const { t, s } = await read('chr1\t0\t10\n', { base: 1, ends: 'inclusive' });
    expect([t!.start[0], t!.end[0]]).toEqual([0, 10]);
    expect(s.issues?.map((i) => i.code)).toEqual(['seq.source.bed.transform']);
    const g = generateBed6(50, names, 3);
    const { file, index } = indexed(g);
    const r = await BedRuntime.open(
      entry('bed', { index: 'i', base: 1, ends: 'inclusive' }),
      context({ f: file, i: index }),
    );
    expect(r.issues?.map((i) => i.code)).toEqual(['seq.source.bed.transform']);
    const first = g.lines.find((l) => l.end > l.begin)!;
    const request = {
      system: names[first.seq],
      start: first.begin,
      end: first.end,
      level: 1,
      shapes: ['ranges' as const],
    };
    const a = await r.read(request, signal);
    const rows = plainRows(a) as { start: number; end: number }[];
    expect(rows.some((x) => x.start === first.begin && x.end === first.end)).toBe(true);
  });

  it('identical lines keep apart by `#2` in file order, the same whole and by region; unnamed ids @start-end', async () => {
    const text = 'chr1\t5\t9\ta\nchr1\t5\t9\ta\nchr1\t5\t9\t.\nchr1\t5\t9\t.\n';
    const { t } = await read(text);
    expect(columnValues(t!.id!)).toEqual(['a@5-9', 'a@5-9#2', '@5-9', '@5-9#2']);
    // by region, the same ids: identical lines meet the same regions, so the suffix is the file's order there too
    const lines = text
      .trim()
      .split('\n')
      .reduce<{ seq: number; begin: number; end: number; start: number; stop: number }[]>((out, l) => {
        const start = out.length ? out[out.length - 1].stop : 0;
        out.push({ seq: 0, begin: 5, end: 9, start, stop: start + l.length + 1 });
        return out;
      }, []);
    const { file, index } = indexed({ text, lines, names: ['chr1'] });
    const r = await BedRuntime.open(entry('bed', { index: 'i' }), context({ f: file, i: index }));
    for (const [start, end] of [
      [0, 10],
      [6, 7],
    ]) {
      const a = await r.read({ system: 'chr1', start, end, level: 1, shapes: ['ranges'] }, signal);
      expect(a.ok && a.tables.ranges?.form === 'list' && columnValues(a.tables.ranges.id!)).toEqual(
        columnValues(t!.id!),
      );
    }
  });

  it('track, browser and # lines read past; lines without a tab split on spaces, said once; CRLF', async () => {
    const { t, s } = await read('browser position chr1\ntrack name=x\n# c\r\nchr1 1 5 a\r\nchr1  6  9  b\r\n');
    expect(Array.from(t!.start)).toEqual([1, 6]);
    expect(s.issues?.map((i) => i.code)).toEqual(['seq.source.bed.whitespace']);
  });

  it('bounds before allocation: bytes and lines of a whole read, pieces over maxRows', async () => {
    const big = generateBed12(300, names, 3).text;
    const small = { maxRows: 100, maxBytes: 1 << 20 };
    // past `maxBytes` refused before a byte is decoded to text
    const decode = vi.spyOn(TextDecoder.prototype, 'decode');
    expect(
      await opened(BedRuntime.open(entry('bed'), context({ f: big }, true, { maxRows: 1e6, maxBytes: 100 }))),
    ).toEqual(['seq.source.too_large']);
    expect(decode).not.toHaveBeenCalled();
    decode.mockRestore();
    // a line of four million blocks past the rows left: refused from its count, its lists never split (they are not
    // even well formed: a parse would have said `.line`)
    const split = vi.spyOn(String.prototype, 'split');
    const huge = `chr1\t0\t100\tx\t0\t+\t0\t100\t0\t4000000\t${'1,'.repeat(1000)}\tnot,a,list`;
    expect(await opened(BedRuntime.open(entry('bed'), context({ f: `${huge}\n` }, true, small)))).toEqual([
      'seq.source.too_large',
    ]);
    expect(split.mock.calls.filter((c) => (c[0] as unknown) === ',').length).toBe(0);
    split.mockRestore();
    expect(await opened(BedRuntime.open(entry('bed'), context({ f: big }, true, small)))).toEqual([
      'seq.source.too_large',
    ]);
    expect(
      await opened(BedRuntime.open(entry('bed'), context({ f: big }, true, { maxRows: 1e6, maxBytes: 100 }))),
    ).toEqual(['seq.source.too_large']);
    // 60 lines of 5 blocks: 300 pieces past 100 rows
    const blocks = Array.from(
      { length: 60 },
      (_, i) => `chr1\t${i * 100}\t${i * 100 + 50}\tx\t0\t+\t0\t0\t0\t5\t10,10,10,10,10,\t0,10,20,30,40,`,
    );
    expect(await opened(BedRuntime.open(entry('bed'), context({ f: blocks.join('\n') }, true, small)))).toEqual([
      'seq.source.too_large',
    ]);
  });

  it('an index the source names is read (EC-186): never `<uri>.tbi` assumed', async () => {
    const g = generateBed12(100, names, 8);
    const { file, index } = indexed(g);
    const files = { f: file, 'f.tbi': index };
    // no `index`: read whole (BGZF by its magic), the `.tbi` beside it never asked
    const whole = await BedRuntime.open(entry('bed'), context(files));
    expect(whole.coverage.levels.length).toBe(2);
    const byRegion = await BedRuntime.open(entry('bed', { index: 'f.tbi' }), context(files));
    expect(byRegion.coverage.levels).toEqual([1]);
    expect(await opened(BedRuntime.open(entry('bed', { index: 'missing.tbi' }), context(files)))).toEqual([
      'seq.source.unreadable',
    ]);
  });
});

describe('bedGraph', () => {
  it('values a position (NaN between), a later interval over an earlier, served as values only', async () => {
    const s = await BedRuntime.open(
      entry('bedgraph'),
      context({ f: 'track type=bedGraph\nchr1\t2\t5\t1.5\nchr1\t4\t6\t3\nchr1\t8\t9\t-2\n' }),
    );
    const a = await s.read({ system: 'chr1', start: 0, end: 20, level: 1, shapes: ['values'] }, signal);
    expect(plainRows(a)).toEqual([0, null, null, 1.5, 1.5, 3, 3, null, null, -2]);
    expect(await s.read({ system: 'chr1', start: 0, end: 20, level: 1, shapes: ['ranges'] }, signal)).toMatchObject({
      ok: false,
      kind: 'not-served',
    });
  });
});

describe('chrom.sizes (EC-188)', () => {
  it('a coverage of systems and lengths, no shapes; a length past int32 refused by name; an index said and ignored', async () => {
    const s = await BedRuntime.open(
      entry('chrom-sizes', { index: 'x' }),
      context({ f: 'chr1\t248956422\nchr2 242193529\naxolotl1\t3000000000\nchr1\t5\nbad\n' }),
    );
    expect(s.coverage).toEqual({
      systems: [
        { name: 'chr1', length: 248956422 },
        { name: 'chr2', length: 242193529 },
      ],
      shapes: [],
      levels: [1],
    });
    expect(s.issues?.map((i) => i.code).sort()).toEqual([
      'seq.source.bed.index',
      'seq.source.bed.line',
      'seq.source.bed.line',
      'seq.source.bed.whitespace',
      'seq.source.too_large',
    ]);
    expect(s.issues?.find((i) => i.code === 'seq.source.too_large')?.message).toContain('axolotl1');
  });
});

describe('chrom.sizes at scale', () => {
  it('20,000 lines read in one pass (a name seen kept in a set, not searched for: 1M in the heavy lane)', async () => {
    const text = Array.from({ length: 20_000 }, (_, i) => `ctg${i}\t${1000 + i}`).join('\n');
    const t0 = performance.now();
    const s = await BedRuntime.open(entry('chrom-sizes'), context({ f: `${text}\nctg5\t9\n` }));
    // generous: quadratic took seconds here (80,000 lines 16.9 s); linear a few tens of ms
    expect(performance.now() - t0).toBeLessThan(1500);
    expect(s.coverage.systems.length).toBe(20_000);
    expect(s.issues?.map((i) => i.message)).toEqual(['A name again (the first kept): 1 (line 20001)']);
  });
});

describe('BEDPE: read whole, indexed by both ends (EC-206, EC-208)', () => {
  const text = generateBedpe(3000, names, 3);

  it('every pair with an end in a region, found from either end, each once; strands on every row', async () => {
    const s = (await BedRuntime.open(entry('bedpe', { index: 'x.tbi' }), context({ f: text }))) as BedpeSource;
    expect(s.issues?.map((i) => i.code)).toEqual(['seq.source.bed.index']);
    expect(s.summaryBytes).toBe(3000 * (2 + 16 + 8 + 2 + 4 + 4));
    const lines = text
      .trim()
      .split('\n')
      .map((l) => l.split('\t'));
    const meets = (sys: string, s0: string, e0: string, system: string, start: number, end: number) =>
      sys === system && Number(s0) < end && Math.max(Number(e0), Number(s0) + 1) > start;
    const rng = new PCG(4);
    for (let k = 0; k < 30; k++) {
      const system = names[k % 3];
      const start = Math.floor(rng.float() * 1_000_000);
      const end = start + Math.floor(rng.float() * 50_000);
      const a = await s.pairs(system, start, end, signal);
      expect(a.ok).toBe(true);
      if (!a.ok) continue;
      const truth = lines
        .map((l, i) => [l, i] as const)
        .filter(([l]) => meets(l[0], l[1], l[2], system, start, end) || meets(l[3], l[4], l[5], system, start, end))
        .map(([, i]) => i);
      expect(Array.from(a.pairs.row!)).toEqual(truth);
      for (let i = 0; i < truth.length; i++) {
        const l = lines[truth[i]];
        expect([a.pairs.aSystem!.get(i), a.pairs.aStart[i], a.pairs.bSystem!.get(i), a.pairs.bEnd[i]]).toEqual([
          l[0],
          Number(l[1]),
          l[3],
          Number(l[5]),
        ]);
        expect(a.pairs.id!.get(i)).toBe(l[6] === '.' ? null : l[6]);
        expect(a.pairs.aStrand![i]).toBe(l[8] === '+' ? 1 : -1);
      }
    }
    // a pair with both ends in the region is found once
    const both = await BedRuntime.open(entry('bedpe'), context({ f: 'chr1\t10\t20\tchr1\t30\t40\tx\t1\t+\t-\n' }));
    const once = await (both as BedpeSource).pairs('chr1', 0, 100, signal);
    expect(once.ok && Array.from(once.pairs.row!)).toEqual([0]);
  });

  it('past `maxPairsBytes` the source is too large, hinting pairix, before the rows are parsed', async () => {
    const r = bedRuntime({ maxPairsBytes: 3000 * 20 });
    try {
      await r.open(entry('bedpe'), context({ f: text }));
      expect.unreachable();
    } catch (e) {
      expect(e).toBeInstanceOf(SourceOpenError);
      const i = (e as SourceOpenError).issues[0];
      expect(i.code).toBe('seq.source.too_large');
      expect(i.hint).toMatch(/pairix/);
    }
  });

  it('12,000 links indexed by both ends (the milestone; 600,000 in the heavy lane)', () => {
    const schema = pairsSchema(10);
    if ('error' in schema) throw new Error(schema.error);
    const b = new PairsBuilder(schema);
    const causes = new BedCauses();
    for (const l of generateBedpe(12_000, names, 9).trim().split('\n')) b.add(l.split('\t'), '', causes);
    const rows = b.finish();
    const index = indexPairs(rows);
    for (let s = 0; s < names.length; s++) {
      const o = index.a.order.subarray(index.a.first[s], index.a.first[s + 1]);
      for (let k = 1; k < o.length; k++) expect(rows.aStart[o[k]] >= rows.aStart[o[k - 1]]).toBe(true);
    }
    expect(pairsIn(rows, index, 0, 0, 2_000_000).length).toBeGreaterThan(0);
  });
});
