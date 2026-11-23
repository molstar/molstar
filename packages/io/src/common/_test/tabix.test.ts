/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { itWithin } from '@molstar/test-harness/timing';
import { bgzfReader, inflateBgzf } from '../bgzf.js';
import { binCount, metaBin, reg2bin, reg2bins } from '../internal/tabix-bins.js';
import {
  readHeaderLines,
  readRegion,
  readTabixIndex,
  regionBytes,
  regionChunks,
  TabixError,
  type TabixIndex,
  type TabixLines,
} from '../tabix.js';
import {
  bgzip,
  type GeneratedFile,
  generateBed,
  generateVcf,
  type IndexOptions,
  memoryRanges,
  tabixIndex,
  tabixIndexRaw,
} from '@molstar/test-harness/genome/bgzf-writer';

const names = ['chr1', 'chr2', 'chrX'];
const bedOptions = (kind: 'tbi' | 'csi'): IndexOptions => ({
  kind,
  format: 0,
  zeroBased: true,
  columns: [1, 2, 3],
  meta: '#',
  skip: 0,
  names,
});
const vcfOptions = (kind: 'tbi' | 'csi'): IndexOptions => ({
  kind,
  format: 2,
  zeroBased: false,
  columns: [1, 2, 0],
  meta: '#',
  skip: 0,
  names,
});

/** A generated file bgzipped in small blocks (lines across block edges) and its index. */
async function fixture(gen: GeneratedFile, options: IndexOptions, blockSize = 2000) {
  const file = bgzip(gen.text, blockSize);
  const index = await readTabixIndex(tabixIndex(gen.lines, file, options));
  const ranges = memoryRanges(file.bytes);
  return { gen, file, index, ranges, reader: bgzfReader(ranges) };
}

const decoder = new TextDecoder();
const linesOf = (r: TabixLines) =>
  Array.from(r.lineStart, (s, i) => ({
    text: decoder.decode(r.bytes.subarray(s, r.lineEnd[i])),
    begin: r.begin[i],
    end: r.end[i],
  }));

/** The whole file's records clipped to a region: the oracle every region read equals. */
function clipped(gen: GeneratedFile, name: string, start: number, end: number) {
  const seq = gen.names.indexOf(name);
  return gen.lines
    .filter((l) => l.seq === seq && l.begin < end && l.end > start)
    .map((l) => ({ text: l.text, begin: l.begin, end: l.end }));
}

describe('bins', () => {
  it("reg2bin by the SAM spec's levels: 0; 1–8; 9–72; 73–584; 585–4680; 4681–37448", () => {
    expect(reg2bin(0, 1, 14, 5)).toBe(4681);
    expect(reg2bin(2 ** 29 - 1, 2 ** 29, 14, 5)).toBe(37448);
    expect(reg2bin(0, 2 ** 14 + 1, 14, 5)).toBe(585);
    expect(reg2bin(0, 2 ** 29, 14, 5)).toBe(0);
    expect(reg2bin(2 ** 26, 2 ** 26 + 2 ** 25, 14, 5)).toBe(2);
    expect(binCount(5)).toBe(37449);
    expect(metaBin(5)).toBe(37450);
  });

  it('reg2bins is every bin whose interval meets the region (brute force), TBI and a CSI scheme', () => {
    const rng = new PCG(5);
    for (const [minShift, depth] of [
      [14, 5],
      [12, 4],
      [14, 6],
    ]) {
      const max = 2 ** (minShift + 3 * depth);
      // every bin's interval
      const intervals: [number, number][] = [];
      for (let l = 0, t = 0; l <= depth; t += 8 ** l, l++) {
        const size = 2 ** (minShift + 3 * (depth - l));
        for (let i = 0; i < 8 ** l; i++) intervals[t + i] = [i * size, (i + 1) * size];
      }
      for (let k = 0; k < (depth > 5 ? 8 : 50); k++) {
        const beg = Math.floor(rng.float() * max);
        const end = Math.min(max, beg + 1 + Math.floor(rng.float() ** 3 * max * 0.01));
        const expected = intervals.flatMap(([a, b], bin) => (a < end && b > beg ? [bin] : []));
        expect(reg2bins(beg, end, minShift, depth)).toEqual(expected);
        // the bin a record is filed in is among the bins of any region it meets
        expect(reg2bins(beg, end, minShift, depth)).toContain(reg2bin(beg, end, minShift, depth));
      }
    }
  });
});

describe('tabix index', () => {
  it('reads a TBI: format, columns, meta, names, scheme', async () => {
    const { index } = await fixture(generateBed(200, names), bedOptions('tbi'));
    expect(index).toMatchObject({
      kind: 'tbi',
      format: 'generic',
      zeroBased: true,
      columns: { seq: 1, begin: 2, end: 3 },
      meta: '#',
      skip: 0,
      names,
      minShift: 14,
      depth: 5,
    });
    const vcf = await fixture(generateVcf(200, names), vcfOptions('csi'));
    expect(vcf.index).toMatchObject({ kind: 'csi', format: 'vcf', zeroBased: false, columns: { end: 0 } });
  });

  it("a chunk starting before the region's least offset begins there (htslib's reg2intervals)", async () => {
    // three lines in one 128 kb bin, one chunk; the region meets only the third's window
    const rows = [
      [10_000, 30_000],
      [50_000, 70_000],
      [110_000, 125_000],
    ];
    const enc = new TextEncoder();
    let at = 0;
    const lines = rows.map(([b, e], i) => {
      const t = `chr1\t${b}\t${e}\tl${i}\n`;
      const line = { seq: 0, begin: b, end: e, start: at, stop: at + t.length, text: t };
      at += t.length;
      return line;
    });
    const file = bgzip(enc.encode(lines.map((l) => l.text).join('')));
    const index = await readTabixIndex(tabixIndex(lines, file, { ...bedOptions('tbi'), names: ['chr1'] }));
    expect(regionChunks(index, 'chr1', 0, 200_000)).toEqual([{ begin: 0, end: file.voffset(at) }]);
    expect(regionChunks(index, 'chr1', 115_000, 116_000)).toEqual([
      { begin: file.voffset(lines[2].start), end: file.voffset(at) },
    ]);
  });

  it('a chunk inside one block costs that block at its bound, never nothing', async () => {
    const { index } = await fixture(generateBed(50, names), bedOptions('tbi'), 65280);
    const chunks = regionChunks(index, 'chr1', 0, 2 ** 29);
    expect(chunks.length).toBe(1);
    expect(Math.floor(chunks[0].begin / 65536)).toBe(Math.floor(chunks[0].end / 65536));
    expect(regionBytes(index, 'chr1', 0, 2 ** 29)).toBe(65536);
  });

  it('regionBytes is the sum of the chunks regionChunks gives', async () => {
    const { index } = await fixture(generateBed(3000, names), bedOptions('tbi'), 1500);
    const rng = new PCG(9);
    for (let k = 0; k < 30; k++) {
      const start = Math.floor(rng.float() * 2_000_000);
      const end = start + Math.floor(rng.float() * 300_000);
      const chunks = regionChunks(index, 'chr2', start, end);
      // each chunk from its first block to its last, the last counted at its bound when read (its size is not in the index)
      const sum = chunks.reduce(
        (n, c) => n + Math.floor(c.end / 65536) - Math.floor(c.begin / 65536) + (c.end % 65536 > 0 ? 65536 : 0),
        0,
      );
      expect(regionBytes(index, 'chr2', start, end)).toBe(sum);
      // merged: sorted, apart
      for (let i = 1; i < chunks.length; i++) expect(chunks[i].begin).toBeGreaterThan(chunks[i - 1].end);
    }
    expect(regionChunks(index, 'chrNone', 0, 1000)).toEqual([]);
  });
});

describe('a region read equals the whole read clipped', () => {
  for (const kind of ['tbi', 'csi'] as const) {
    it(`BED by ${kind}: random regions, lines across block edges`, async () => {
      const f = await fixture(generateBed(4000, names, 3), bedOptions(kind), 1700);
      const rng = new PCG(13);
      for (let k = 0; k < 60; k++) {
        const name = names[k % names.length];
        const start = Math.floor(rng.float() * 2_100_000);
        const end = start + 1 + Math.floor(rng.float() ** 2 * 400_000);
        const got = await readRegion(f.reader, f.index, name, start, end);
        expect(linesOf(got)).toEqual(clipped(f.gen, name, start, end));
      }
      // the whole sequence
      const all = await readRegion(f.reader, f.index, 'chrX', 0, 2 ** 29);
      expect(linesOf(all)).toEqual(clipped(f.gen, 'chrX', 0, 2 ** 29));
    });

    it(`VCF by ${kind}: extents by tbx_parse1 (REF's length, INFO/END on any line, not SVEND)`, async () => {
      const f = await fixture(generateVcf(3000, names, 4), vcfOptions(kind), 1500);
      const rng = new PCG(17);
      for (let k = 0; k < 60; k++) {
        const name = names[k % names.length];
        const start = Math.floor(rng.float() * 1_000_000);
        const end = start + 1 + Math.floor(rng.float() ** 2 * 200_000);
        expect(linesOf(await readRegion(f.reader, f.index, name, start, end))).toEqual(
          clipped(f.gen, name, start, end),
        );
      }
      const all = linesOf(await readRegion(f.reader, f.index, 'chr1', 0, 2 ** 29));
      expect(all.some((l) => /^chr1\t\d+\tid\d+\tA\tT\t50\tPASS\tEND=/.test(l.text) && l.end - l.begin === 101)).toBe(
        true,
      );
      expect(all.some((l) => /SVEND=/.test(l.text) && l.end - l.begin === 1)).toBe(true);
    });
  }

  it('CSI past 512 Mb (a TBI cannot index it): positions near 600 Mb read', async () => {
    const gen = generateBed(500, ['chrBig'], 21, 3000, 600_000_000);
    const options: IndexOptions = { ...bedOptions('csi'), names: ['chrBig'], depth: 6 };
    const f = await fixture(gen, options);
    const start = 600_500_000;
    const end = 601_000_000;
    const got = linesOf(await readRegion(f.reader, f.index, 'chrBig', start, end));
    expect(got.length).toBeGreaterThan(10);
    expect(got).toEqual(clipped(gen, 'chrBig', start, end));
  });

  it('reads the header lines, and only them', async () => {
    const f = await fixture(generateVcf(100, names), vcfOptions('tbi'), 50);
    const header = linesOf(await readHeaderLines(f.reader, f.index));
    expect(header.map((l) => l.text)).toEqual(f.gen.header);
    expect(header.every((l) => l.begin === -1 && l.end === -1)).toBe(true);
  });
});

describe('bounds and malformed files', () => {
  const tooLarge = { name: 'TabixError', kind: 'too-large' };
  const malformed = { name: 'TabixError', kind: 'malformed' };

  it('a region past maxRegionBytes is too-large before a byte is read', async () => {
    const f = await fixture(generateBed(4000, names), bedOptions('tbi'), 1000);
    const bytes = regionBytes(f.index, 'chr1', 0, 2 ** 29);
    expect(bytes).toBeGreaterThan(10_000);
    await expect(
      readRegion(f.reader, f.index, 'chr1', 0, 2 ** 29, { maxRegionBytes: bytes - 1 }),
    ).rejects.toMatchObject(tooLarge);
    expect(f.ranges.reads).toBe(0);
  });

  it('a line past maxLineBytes is malformed; an index past maxIndexBytes too-large', async () => {
    const f = await fixture(generateBed(400, names), bedOptions('tbi'));
    await expect(readRegion(f.reader, f.index, 'chr1', 0, 2 ** 29, { maxLineBytes: 10 })).rejects.toMatchObject(
      malformed,
    );
    const raw = tabixIndexRaw(f.gen.lines, f.file, bedOptions('tbi'));
    await expect(readTabixIndex(bgzip(raw).bytes, { maxIndexBytes: raw.length - 1 })).rejects.toMatchObject(tooLarge);
    await expect(readTabixIndex(raw, { maxIndexBytes: raw.length - 1 })).rejects.toMatchObject(tooLarge);
  });

  it('counts past the index bytes are refused before they allocate', async () => {
    const f = await fixture(generateBed(50, names), bedOptions('tbi'));
    const raw = tabixIndexRaw(f.gen.lines, f.file, bedOptions('tbi'));
    const at = (o: number, v: number) => {
      const c = raw.slice();
      new DataView(c.buffer).setInt32(o, v, true);
      return c;
    };
    const lnm = new DataView(raw.buffer).getInt32(32, true);
    const firstRef = 36 + lnm;
    for (const bad of [
      at(4, 2 ** 30), // n_ref
      at(4, -1),
      at(4, 2), // names for 3 references
      at(8, 7), // format
      at(12, 0), // a column 0
      at(32, 2 ** 30), // l_nm
      at(firstRef, 2 ** 30), // n_bin
      at(firstRef, 40000),
      at(firstRef + 4, 50000), // a bin past the scheme
      at(firstRef + 8, 2 ** 29), // n_chunk
      raw.subarray(0, raw.length - 30),
      new TextEncoder().encode('TBX\x01'),
    ]) {
      await expect(readTabixIndex(bad)).rejects.toMatchObject(malformed);
    }
    // a chunk ending before it begins
    const swapped = raw.slice();
    const v = new DataView(swapped.buffer);
    const c0 = firstRef + 12;
    const b = v.getBigUint64(c0, true);
    v.setBigUint64(c0, v.getBigUint64(c0 + 8, true), true);
    v.setBigUint64(c0 + 8, b, true);
    await expect(readTabixIndex(swapped)).rejects.toMatchObject(malformed);
  });

  it('an index for another file gives malformed regions, never rows in the wrong place', async () => {
    const a = await fixture(generateBed(2000, names, 1), bedOptions('tbi'), 1500);
    const b = await fixture(generateBed(2000, names, 99), bedOptions('tbi'), 1500);
    const rng = new PCG(23);
    let refused = 0;
    for (let k = 0; k < 30; k++) {
      const name = names[k % 3];
      const start = Math.floor(rng.float() * 2_000_000);
      const end = start + 1 + Math.floor(rng.float() * 300_000);
      try {
        const got = linesOf(await readRegion(a.reader, b.index, name, start, end));
        // whatever is answered is a's lines of this region, with their right extents
        const truth = clipped(a.gen, name, start, end);
        for (const l of got) expect(truth).toContainEqual(l);
      } catch (e) {
        expect((e as Error).name).toMatch(/TabixError|BgzfError/);
        refused++;
      }
    }
    expect(refused).toBeGreaterThan(0);
  });

  it('fuzzed and truncated indexes never throw past a TabixError', async () => {
    const f = await fixture(generateVcf(300, names), vcfOptions('tbi'));
    const raw = tabixIndexRaw(f.gen.lines, f.file, vcfOptions('tbi'));
    const csi = tabixIndexRaw(f.gen.lines, f.file, vcfOptions('csi'));
    const rng = new PCG(29);
    for (let k = 0; k < 300; k++) {
      const src = k % 2 ? raw : csi;
      const c = src.slice(0, k % 5 === 0 ? Math.floor(rng.float() * src.length) : src.length);
      for (let j = 0; j < 1 + (k % 3); j++) c[Math.floor(rng.float() * c.length)] = Math.floor(rng.float() * 256);
      let index: TabixIndex | undefined;
      try {
        index = await readTabixIndex(c);
      } catch (e) {
        expect((e as Error).name, `${e}`).toBe('TabixError');
      }
      if (!index) continue;
      for (const name of index.names.slice(0, 2)) {
        try {
          await readRegion(f.reader, index, name, 0, 2 ** 29);
        } catch (e) {
          expect((e as Error).name, `${e}`).toMatch(/^(TabixError|BgzfError)$/);
        }
      }
    }
    expect(TabixError.name).toBe('TabixError');
  });
});

describe('header lines with no line indexed', () => {
  it('scans each block once: allocation bounded by the header, not blocks × header', async () => {
    const header = Array.from({ length: 4000 }, (_, i) => `##meta line ${i} ${'x'.repeat(40)}`);
    const text = new TextEncoder().encode(header.map((h) => `${h}\n`).join(''));
    const file = bgzip(text, 500);
    const index = await readTabixIndex(tabixIndex([], file, vcfOptions('tbi')));
    const reader = bgzfReader(memoryRanges(file.bytes), { cacheBlocks: 1 });
    let allocated = 0;
    const Native = globalThis.Uint8Array;
    vi.stubGlobal(
      'Uint8Array',
      class extends Native {
        constructor(...args: ConstructorParameters<typeof Native>) {
          super(...args);
          allocated += this.byteLength;
        }
      },
    );
    let got: TabixLines;
    try {
      got = await readHeaderLines(reader, index);
    } finally {
      vi.unstubAllGlobals();
    }
    expect(got.lineStart.length).toBe(header.length);
    // the header is about 200 kB in about 400 blocks: joining all after each block is tens of MB
    expect(allocated).toBeLessThan(8 * text.length);
  });
});

describe('performance', () => {
  // A region read, quick here, bounded in the heavy lane
  const gen = generateVcf(30_000, ['chr1'], 31);
  const file = bgzip(gen.text);
  const indexBytes = tabixIndex(gen.lines, file, { ...vcfOptions('tbi'), names: ['chr1'] });
  itWithin('a region of about 10,000 VCF lines read by tabix', async (within) => {
    const index = await readTabixIndex(indexBytes);
    const reader = bgzfReader(memoryRanges(file.bytes));
    const got = await within(100, () => readRegion(reader, index, 'chr1', 0, 333_000));
    expect(got.lineStart.length).toBeGreaterThan(9000);
    expect(got.lineStart.length).toBe(clipped(gen, 'chr1', 0, 333_000).length);
  });

  it('the whole file inflated equals the generated text', async () => {
    expect(await inflateBgzf(file.bytes)).toEqual(gen.text);
  });
});
