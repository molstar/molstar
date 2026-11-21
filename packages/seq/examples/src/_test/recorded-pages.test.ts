/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { describe, expect, it } from 'vitest';
import { nodeFixtureAssets, PackageFolder } from '../node-assets.js';
import {
  encodePages,
  PageNotRecorded,
  PageSize,
  PagesSuffix,
  pagesOf,
  readRecording,
  recordedAssets,
  recordedRanges,
  recordingRanges,
  withRecordedPages,
} from '../recorded-pages.js';

// Real files as recorded pages: a recording reads back byte for byte (the hg38 cytobands' checksum is the download's,
// in `sources.json`), a page not recorded fails naming it, the recorder notes the pages a reader asks, and a runtime
// reads its file from the recording through the asset port.

const signal = new AbortController().signal;
const Cytobands = 'fixtures/genomes/cytoBandIdeo.txt.gz';

describe('recorded pages', () => {
  it('the cytobands read back as downloaded: the checksum sources.json gives', async () => {
    const assets = recordedAssets(nodeFixtureAssets().reader, [Cytobands]);
    const bytes = await assets.readBytes(Cytobands, signal);
    const rows = JSON.parse(readFileSync(join(PackageFolder, 'fixtures/genomes/sources.json'), 'utf8'));
    expect(createHash('sha256').update(bytes).digest('hex')).toBe(rows[0].sha256Download);
    expect(bytes.length).toBe(rows[0].bytesDownloaded);
  });

  it('a read of a page not recorded fails naming the page and the file', async () => {
    const file = new Uint8Array(3 * PageSize + 10).map((_, i) => i % 251);
    const r = recordedRanges(encodePages('big.bin', file, [0, 2]));
    const [a, b] = await r.read([
      { offset: 5, length: 10 },
      { offset: 2 * PageSize - 3, length: 0 },
    ]);
    expect(Array.from(a)).toEqual(Array.from(file.subarray(5, 15)));
    expect(b.length).toBe(0);
    const across = (await r.read([{ offset: 2 * PageSize + 1, length: 20 }]))[0];
    expect(Array.from(across)).toEqual(Array.from(file.subarray(2 * PageSize + 1, 2 * PageSize + 21)));
    await expect(r.read([{ offset: PageSize - 2, length: 4 }])).rejects.toThrow(PageNotRecorded);
    await expect(r.read([{ offset: 3 * PageSize, length: 4 }])).rejects.toThrow(/Page 3 of big\.bin/);
  });

  it('the recorder notes the pages a reader asks; they are what a recording needs', async () => {
    const file = new Uint8Array(5 * PageSize);
    const rec = recordingRanges(file);
    await rec.read([
      { offset: PageSize + 5, length: PageSize },
      { offset: 4 * PageSize, length: 1 },
    ]);
    expect([...rec.pages].sort()).toEqual([1, 2, 4]);
    expect(pagesOf([{ offset: 0, length: PageSize }], file.length)).toEqual([0]);
    const replay = recordedRanges(encodePages('f', file, rec.pages));
    await expect(replay.read([{ offset: PageSize + 5, length: PageSize }])).resolves.toHaveLength(1);
  });

  it('a runtime reads its file from the recording: hg38’s cytobands as BED', async () => {
    const runtime = withRecordedPages(BedRuntime, [Cytobands]);
    const s = await runtime.open(
      { address: '/data/bands', format: 'bed', uri: Cytobands, fields: { stain: '4' } },
      { assets: nodeFixtureAssets().reader, signal, limits: { maxRows: 1e6, maxBytes: 64 << 20 } },
    );
    expect(s.issues).toEqual([]);
    expect(s.coverage.systems.map((x) => x.name)).toContain('chr17');
  });

  it('a page cut short is refused naming it and the file, never read as zeros', () => {
    const file = new Uint8Array(PageSize + 100).fill(7);
    const r = encodePages('f.bin', file, 'all');
    const cut = { ...r, pages: { ...r.pages, 1: r.pages[1].slice(0, 100) } };
    expect(() => readRecording(JSON.stringify(cut))).toThrow(/Page 1 of f\.bin is not its 100 bytes/);
    const short = { ...r, pages: { ...r.pages, 0: r.pages[1] } };
    expect(() => readRecording(JSON.stringify(short))).toThrow(/Page 0 of f\.bin/);
    expect(() => readRecording(JSON.stringify({ ...r, pages: { 5: r.pages[1] } }))).toThrow(/Page 5 of f\.bin/);
    expect(readRecording(JSON.stringify(r)).size).toBe(PageSize + 100);
  });

  it('only declared files are asked for a recording; a bad recording or an abort fails, never falls back', async () => {
    const asked: string[] = [];
    const files: Record<string, string> = {
      'plain.txt': 'plain',
      'f.bin': 'the file itself, never read',
      [`f.bin${PagesSuffix}`]: '{"format":"seqx-pages/1","file":"f.bin","size":3,"pageSize":65536,"pages":{"0":"AA"}}',
    };
    const inner = {
      readText: async (uri: string, s: AbortSignal) => {
        asked.push(uri);
        s.throwIfAborted();
        if (!(uri in files)) throw new Error(`no ${uri}`);
        return files[uri];
      },
      readBytes: async (uri: string) => {
        asked.push(uri);
        return new TextEncoder().encode(files[uri]);
      },
    };
    const assets = recordedAssets(inner, ['f.bin', 'missing.bin']);
    expect(await assets.readText('plain.txt', signal)).toBe('plain');
    expect(asked).toEqual(['plain.txt']);
    await expect(assets.readBytes('f.bin', signal)).rejects.toThrow(/Page 0 of f\.bin/);
    await expect(assets.readBytes('missing.bin', signal)).rejects.toThrow(/no missing\.bin\.pages\.json/);
    const stop = new AbortController();
    stop.abort();
    await expect(assets.readBytes('missing.bin', stop.signal)).rejects.toMatchObject({ name: 'AbortError' });
    expect(asked).not.toContain('f.bin');
    expect(asked).not.toContain('missing.bin');
  });
});
