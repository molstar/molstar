/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { existsSync, mkdtempSync, readdirSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { describe, expect, it } from 'vitest';
import { cutLines, isPlainName, parseCut, recordPages } from '../record-pages.ts';

// The recording script never writes outside the fixtures folder: a dataset or a file name is one plain path segment,
// refused before anything is read, fetched or written. A cut keeps a region's lines and the header.

describe('record-pages', () => {
  it('refuses a dataset or a name that leaves the fixtures folder, writing nothing', async () => {
    const root = mkdtempSync(join(tmpdir(), 'pages-'));
    const fixtures = join(root, 'fixtures');
    const from = join(root, 'download.bin');
    writeFileSync(from, 'abc');
    const job = {
      url: 'https://example.org/x.bin',
      entry: 'e',
      licence: 'CC0-1.0',
      attribution: 'a',
      pages: 'all' as const,
      from,
    };
    for (const dataset of ['../out', '..', 'a/b', '/abs', '.hidden', 'a\\b'])
      await expect(recordPages({ ...job, dataset }, fixtures)).rejects.toThrow(/not a plain folder name/);
    for (const name of ['../x', 'a/b', '..'])
      await expect(recordPages({ ...job, dataset: 'ok', name }, fixtures)).rejects.toThrow(/not a plain file name/);
    expect(existsSync(fixtures)).toBe(false);
    expect(readdirSync(root)).toEqual(['download.bin']);
    // a plain one is written in its dataset's folder
    expect(isPlainName('genomes')).toBe(true);
    await recordPages({ ...job, dataset: 'ok' }, fixtures);
    expect(readdirSync(join(fixtures, 'ok')).sort()).toEqual(['sources.json', 'x.bin.pages.json']);
  });

  it('a cut keeps the header (its sequence’s `##sequence-region` only) and the lines meeting the region', async () => {
    const lines = [
      '##gff-version 3',
      '##sequence-region 1 1 1000',
      '##sequence-region 17 1 83257441',
      '1\t.\tgene\t5\t10\t.\t+\t.\tID=a',
      '17\t.\tgene\t90\t99\t.\t+\t.\tID=before',
      '17\t.\tgene\t90\t100\t.\t+\t.\tID=edge',
      '###',
      '17\t.\texon\t150\t160\t.\t+\t.\tParent=in',
      '17\t.\tgene\t201\t300\t.\t+\t.\tID=after',
      '##FASTA',
      '17\t.\tgene\t150\t160\t.\t+\t.\tID=fasta',
    ];
    const cut = parseCut('17:100-200');
    expect(cut).toEqual({ seq: '17', start: 100, end: 200 });
    expect(parseCut('chr17:7,668,402-7,687,550')).toEqual({ seq: 'chr17', start: 7668402, end: 7687550 });
    expect(() => parseCut('17:200-100')).toThrow();
    expect(await cutLines(lines, cut)).toEqual([lines[0], lines[2], lines[5], lines[7]]);
  });
});
