/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { columnValues } from '@molstar/core/repr/tree';
import type { AssetRangeReader } from '@molstar/core/util/assets';
import { describe, expect, it } from 'vitest';
import { sourceSuite } from '../../data/_test/source-suite.js';
import { type AssetReader, DefaultOpenLimits, type OpenContext } from '../../data/protocol.js';
import { VcfRuntime } from '../vcf.js';

// The recorded real slice (_spec/README.md "Data", VCF): NCBI ClinVar's GRCh38 `clinvar_20240107.vcf.gz` (BGZF, its
// header's page, TP53's chunks and the EOF block's page) and its `.tbi` (whole), recorded as pages (examples'
// `fixtures/genomes/`); read by region as a host reads it, a page the recording did not make failing its read.

const signal = new AbortController().signal;
const Fixtures = '../../../../examples/fixtures/genomes/';
const File = 'clinvar_20240107.vcf.gz';

interface Recording {
  readonly size: number;
  readonly pageSize: number;
  readonly pages: Record<string, string>;
}

const recording = (name: string): Recording =>
  JSON.parse(readFileSync(fileURLToPath(new URL(`${Fixtures}${name}.pages.json`, import.meta.url)), 'utf8'));

/** A range reader over a recording: a read of a page it did not make throws, naming the page. */
function recordedRanges(r: Recording): AssetRangeReader {
  const page = (p: number) => {
    const b64 = r.pages[p];
    if (b64 === undefined) throw new Error(`page ${p} not recorded`);
    return new Uint8Array(Buffer.from(b64, 'base64'));
  };
  return {
    size: r.size,
    read: async (ranges) =>
      ranges.map(({ offset, length }) => {
        const end = Math.min(r.size, offset + length);
        const out = new Uint8Array(Math.max(0, end - offset));
        for (let p = Math.floor(offset / r.pageSize); p * r.pageSize < end; p++) {
          const bytes = page(p);
          const from = Math.max(offset, p * r.pageSize);
          const to = Math.min(end, (p + 1) * r.pageSize);
          out.set(bytes.subarray(from - p * r.pageSize, to - p * r.pageSize), from - offset);
        }
        return out;
      }),
    dispose() {},
  };
}

function wholeBytes(r: Recording): Uint8Array {
  const out = new Uint8Array(r.size);
  for (const [k, b64] of Object.entries(r.pages)) out.set(Buffer.from(b64, 'base64'), Number(k) * r.pageSize);
  return out;
}

const file = recording(File);
const index = wholeBytes(recording(`${File}.tbi`));
const assets: AssetReader = {
  readText: async () => {
    throw new Error('not read as text');
  },
  readBytes: async (uri) => {
    if (uri === 'clinvar.tbi') return index;
    throw new Error(`${uri} is not read whole`);
  },
  openRange: async () => recordedRanges(file),
};
const context: OpenContext = { assets, signal, limits: DefaultOpenLimits };
const entry = { address: '/data/clinvar', format: 'vcf', uri: 'clinvar.vcf.gz', index: 'clinvar.tbi' };

/** TP53 in GRCh38 (0-based, end exclusive): chr17:7,668,402-7,687,550, ClinVar's `17`. */
const Tp53 = { system: '17', start: 7_668_401, end: 7_687_550 };

describe('the recorded ClinVar slice (TP53, by region)', () => {
  it('passes the source suite', async () => {
    await sourceSuite(VcfRuntime, entry, context);
  });

  it('its header: no contig lines (the names the index’s: `17`, read onto chr17 by the names rule), ClinVar’s INFO keys typed', async () => {
    const s = await VcfRuntime.open(entry, context);
    expect(s.issues).toEqual([]);
    expect(s.coverage.levels).toEqual([1]);
    expect(s.coverage.systems.map((x) => x.name)).toContain('17');
    expect(s.coverage.systems.every((x) => x.length === undefined)).toBe(true);
    const fields = Object.fromEntries((s.coverage.shapes[0].fields ?? []).map((f) => [f.name, f.type]));
    expect(fields).toMatchObject({ ALLELEID: 'int32', AF_EXAC: 'float32', CLNSIG: 'category', GENEINFO: 'category' });
    s.dispose();
  });

  it('TP53’s variants: a row an allele, R175H among them (EC-196), every row TP53’s or a neighbour’s', async () => {
    const s = await VcfRuntime.open(entry, context);
    const a = await s.read({ ...Tp53, level: 1, shapes: ['ranges'] }, signal);
    expect(a.ok).toBe(true);
    if (!a.ok || a.tables.ranges?.form !== 'list') return;
    const r = a.tables.ranges;
    const col = (name: string) => columnValues(r.fields!.find((f) => f.name === name)!.column);
    const ids = columnValues(r.id!) as string[];
    expect(ids.length).toBeGreaterThan(2000);
    // TP53 c.524G>A, p.Arg175His: chr17:7,675,088 C>T on the forward strand, ClinVar's variation 12374
    const k = ids.indexOf('7675088:C>T');
    expect(k).toBeGreaterThanOrEqual(0);
    expect([r.start[k], r.end[k], col('name')[k], col('category')[k]]).toEqual([7675087, 7675088, '12374', 'SNV']);
    expect(String(col('GENEINFO')[k])).toMatch(/^TP53:7157/);
    expect(new Set(col('category'))).toEqual(new Set(['SNV', 'deletion', 'insertion', 'indel', 'MNV']));
    // every row meets the region, in order of start
    for (let i = 0; i < r.start.length; i++) {
      expect(r.start[i] < Tp53.end && r.end[i] > Tp53.start).toBe(true);
      if (i) expect(r.start[i]).toBeGreaterThanOrEqual(r.start[i - 1]);
    }
    s.dispose();
  });
});
