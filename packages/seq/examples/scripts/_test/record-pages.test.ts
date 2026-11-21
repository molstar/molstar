/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { existsSync, mkdtempSync, readdirSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { describe, expect, it } from 'vitest';
import { isPlainName, recordPages } from '../record-pages.ts';

// The recording script never writes outside the fixtures folder: a dataset or a file name is one plain path segment,
// refused before anything is read, fetched or written.

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
});
