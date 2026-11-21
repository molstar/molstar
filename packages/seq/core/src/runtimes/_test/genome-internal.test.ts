/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { gzipSync } from 'node:zlib';
import { describe, expect, it } from 'vitest';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { type AssetReader, DefaultOpenLimits } from '../../data/protocol.js';
import { SourceOpenError } from '../../data/seq-runtime.js';
import { openIndexedLines, wholeLines } from '../internal/lines.js';
import { decodeText, textBytes } from '../internal/text.js';
import {
  bgzip,
  generateBed,
  type GeneratedFile,
  type IndexOptions,
  join,
  memoryRanges,
  tabixIndex,
} from './genome-writer.js';

const names = ['chr1', 'chr2'];
const bed: IndexOptions = { kind: 'tbi', format: 0, zeroBased: true, columns: [1, 2, 3], meta: '#', skip: 0, names };
const signal = new AbortController().signal;
const decoder = new TextDecoder();

/** An asset port over files in memory, with or without range reads, counting what it is asked. */
function assets(files: Record<string, Uint8Array>, ranges = true) {
  const asked: string[] = [];
  const reader: AssetReader = {
    readText: async (uri) => decoder.decode(files[uri]),
    readBytes: async (uri) => {
      asked.push(`bytes ${uri}`);
      if (!files[uri]) throw new Error('no such file');
      return files[uri];
    },
    ...(ranges
      ? {
          openRange: async (uri: string) => {
            asked.push(`range ${uri}`);
            if (!files[uri]) throw new Error('no such file');
            return memoryRanges(files[uri]);
          },
        }
      : {}),
  };
  return { reader, asked };
}

/** The data lines of a whole read meeting a region, by the generator's extents: the oracle. */
function clipped(gen: GeneratedFile, name: string, start: number, end: number) {
  const seq = gen.names.indexOf(name);
  return gen.lines.filter((l) => l.seq === seq && l.begin < end && l.end > start).map((l) => l.text);
}

describe('text: bytes to text by magic, bounded', () => {
  const text = new TextEncoder().encode('#h\nchr1\t1\t2\n'.repeat(500));

  it('plain, gzip and BGZF read alike; never by extension', async () => {
    for (const [bytes, compression] of [
      [text, 'none'],
      [new Uint8Array(gzipSync(text)), 'gzip'],
      [bgzip(text, 1000).bytes, 'bgzf'],
    ] as const) {
      const r = await textBytes(bytes, 1 << 20);
      expect(r).toMatchObject({ ok: true, compression });
      if (r.ok) expect(r.bytes).toEqual(text);
    }
  });

  it('bounded before and after inflating; a damaged gzip malformed, not logged', async () => {
    expect(await textBytes(text, text.length - 1)).toMatchObject({ ok: false, kind: 'too-large' });
    expect(await textBytes(new Uint8Array(gzipSync(text)), text.length - 1)).toMatchObject({
      ok: false,
      kind: 'too-large',
    });
    expect(await textBytes(bgzip(text, 1000).bytes, text.length - 1)).toMatchObject({ ok: false, kind: 'too-large' });
    const gz = new Uint8Array(gzipSync(text));
    gz[gz.length - 8] ^= 1;
    expect(await textBytes(gz, 1 << 20)).toMatchObject({ ok: false, kind: 'malformed' });
    const cut = bgzip(text, 1000).bytes;
    expect(await textBytes(cut.subarray(0, cut.length - 40), 1 << 20)).toMatchObject({ ok: false, kind: 'malformed' });
  });

  it('decodes UTF-8 and tells a replacement character', () => {
    expect(decodeText(new TextEncoder().encode('ačb'))).toEqual({ text: 'ačb', replaced: false });
    expect(decodeText(new Uint8Array([97, 0xff, 98])).replaced).toBe(true);
  });
});

describe('lines read whole', () => {
  it('header lines at the top (meta, skip), blank lines left out, numbers 1-based, CRLF', () => {
    const bytes = new TextEncoder().encode('skipped\n#a\n#b\r\nx\t1\n\n#c\ny\t2');
    const r = wholeLines(bytes, { meta: '#', skip: 1 });
    expect(r.ok).toBe(true);
    if (!r.ok) return;
    const text = (s: { bytes: Uint8Array; lineStart: Uint32Array; lineEnd: Uint32Array }) =>
      Array.from(s.lineStart, (a, i) => decoder.decode(s.bytes.subarray(a, s.lineEnd[i])));
    expect(text(r.value.header)).toEqual(['skipped', '#a', '#b']);
    expect(text(r.value.lines)).toEqual(['x\t1', '#c', 'y\t2']);
    expect(Array.from(r.value.lines.lineNumber)).toEqual([4, 6, 7]);
  });

  it('refuses a line over its bound and lines over maxLines before allocating them', () => {
    const bytes = new TextEncoder().encode('#h\na\nb\nc\n');
    expect(wholeLines(bytes, { meta: '#', maxLines: 2 })).toMatchObject({ ok: false, kind: 'too-large' });
    expect(wholeLines(bytes, { meta: '#', maxLines: 3 }).ok).toBe(true);
    expect(wholeLines(new TextEncoder().encode('#h\nabcdef\n'), { meta: '#', maxLineBytes: 3 })).toMatchObject({
      ok: false,
      kind: 'malformed',
    });
  });
});

describe('lines by region through tabix', () => {
  const gen = generateBed(3000, names, 5);
  const file = bgzip(gen.text, 1500);
  const index = tabixIndex(gen.lines, file, bed);

  for (const ranges of [true, false]) {
    it(`a region read equals the whole read clipped (${ranges ? 'range reads' : 'the file read whole'})`, async () => {
      const a = assets({ 'f.bed.gz': file.bytes, 'f.bed.gz.tbi': index }, ranges);
      const lines = await openIndexedLines('f.bed.gz', 'f.bed.gz.tbi', {
        assets: a.reader,
        signal,
        limits: DefaultOpenLimits,
        address: 'genes',
      });
      expect(lines.eof).toBe(true);
      expect(a.asked).toEqual(
        ranges ? ['range f.bed.gz', 'bytes f.bed.gz.tbi'] : ['bytes f.bed.gz', 'bytes f.bed.gz.tbi'],
      );
      const rng = new PCG(3);
      for (let k = 0; k < 20; k++) {
        const name = names[k % 2];
        const start = Math.floor(rng.float() * 2_000_000);
        const end = start + 1 + Math.floor(rng.float() * 200_000);
        const got = await lines.region(name, start, end, signal);
        expect(got.ok).toBe(true);
        if (!got.ok) continue;
        const l = got.lines;
        expect(Array.from(l.lineStart, (s, i) => decoder.decode(l.bytes.subarray(s, l.lineEnd[i])))).toEqual(
          clipped(gen, name, start, end),
        );
        expect(lines.cost(name, start, end)).toBeGreaterThanOrEqual(0);
      }
      const header = await lines.header(signal);
      expect(header.ok && header.lines.lineStart.length).toBe(gen.header.length);
      lines.dispose();
    });
  }

  it('a region over its bound answers too-large; an aborted read rejects with an AbortError', async () => {
    const a = assets({ f: file.bytes, i: index });
    const lines = await openIndexedLines('f', 'i', {
      assets: a.reader,
      signal,
      limits: DefaultOpenLimits,
      address: 'genes',
      tabix: { maxRegionBytes: 1000 },
    });
    expect(await lines.region('chr1', 0, 2 ** 29, signal)).toMatchObject({ ok: false, kind: 'too-large' });
    const ac = new AbortController();
    ac.abort();
    await expect(lines.region('chr1', 0, 1000, ac.signal)).rejects.toMatchObject({ name: 'AbortError' });
  });

  it('a file without its EOF block is read, told; an index for another file answers malformed, never rows', async () => {
    const other = generateBed(3000, names, 77);
    const otherIndex = tabixIndex(other.lines, bgzip(other.text, 1500), bed);
    const a = assets({ f: bgzip(gen.text, 1500, false).bytes, i: otherIndex });
    const lines = await openIndexedLines('f', 'i', {
      assets: a.reader,
      signal,
      limits: DefaultOpenLimits,
      address: 'g',
    });
    expect(lines.eof).toBe(false);
    const kinds = new Set<string>();
    for (let k = 0; k < 10; k++) {
      const got = await lines.region('chr2', k * 200_000, k * 200_000 + 150_000, signal);
      if (got.ok) {
        const truth = clipped(gen, 'chr2', k * 200_000, k * 200_000 + 150_000);
        const l = got.lines;
        for (let i = 0; i < l.lineStart.length; i++) {
          expect(truth).toContain(decoder.decode(l.bytes.subarray(l.lineStart[i], l.lineEnd[i])));
        }
      } else kinds.add(got.kind);
    }
    expect(kinds.has('malformed')).toBe(true);
  });

  it('a file or index that cannot be read is a SourceOpenError with its code', async () => {
    const open = (files: Record<string, Uint8Array>) =>
      openIndexedLines('f', 'i', { assets: assets(files).reader, signal, limits: DefaultOpenLimits, address: 'g' });
    const code = async (files: Record<string, Uint8Array>) => {
      try {
        await open(files);
        return 'opened';
      } catch (e) {
        expect(e).toBeInstanceOf(SourceOpenError);
        return (e as SourceOpenError).issues[0].code;
      }
    };
    expect(await code({ i: index })).toBe('seq.source.unreadable');
    expect(await code({ f: file.bytes })).toBe('seq.source.unreadable');
    expect(await code({ f: file.bytes, i: join([new TextEncoder().encode('TBX\x01'), new Uint8Array(40)]) })).toBe(
      'seq.source.malformed',
    );
    const whole = assets({ f: file.bytes, i: index }, false);
    await expect(
      openIndexedLines('f', 'i', {
        assets: whole.reader,
        signal,
        limits: { ...DefaultOpenLimits, maxBytes: 100 },
        address: 'g',
      }),
    ).rejects.toMatchObject({ issues: [{ code: 'seq.source.too_large' }] });
  });
});
