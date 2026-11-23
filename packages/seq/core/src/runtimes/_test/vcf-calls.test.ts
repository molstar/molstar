/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The genotype card says the call (_spec/README.md "Data", VCF, The matrix; "Render", Hover texts by region): the
// sample's call as the VCF wrote it, with a word for what the call is, never a per-allele class written as a call.

import { describe, expect, it } from 'vitest';
import { DefaultOpenLimits, type MatrixItems, type RegionAnswer, type RegionRequest } from '../../data/protocol.js';
import { callWords, matrixCellHover } from '../../render/region-hover.js';
import { VcfRuntime } from '../vcf.js';
import { memoryAssets } from './bed-fixtures.js';

const signal = new AbortController().signal;

/** One call a sample, over two sites: a biallelic one and one with two ALT alleles. */
const Calls = ['0/0', '0/1', '1/1', '0|1', '1/2', '0', '1', './.', './1'];
const text = [
  '##fileformat=VCFv4.2',
  '##contig=<ID=chr1,length=1000>',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${Calls.map((_, i) => `S${i}`).join('\t')}`,
  `chr1\t10\t.\tA\tC,G\t.\t.\t.\tGT\t${Calls.join('\t')}`,
].join('\n');

async function matrix() {
  const s = await VcfRuntime.open(
    { address: '/data/s', format: 'vcf', uri: 'f' },
    { assets: memoryAssets({ f: `${text}\n` }), signal, limits: DefaultOpenLimits },
  );
  const r: RegionRequest = { system: 'chr1', start: 0, end: 100, level: 1, shapes: ['matrix'] };
  const a = (await s.read(r, signal)) as RegionAnswer & { ok: true };
  return a.tables.matrix as MatrixItems;
}

/** The card of sample `s` at the table's first item. */
function cardOf(t: MatrixItems, s: number) {
  const lay = { states: t.states, names: Calls.map((_, i) => `S${i}`), system: 'chr1', level: 1, tables: [t] };
  return matrixCellHover(lay as never, t.start[0], s, (x) => x);
}

describe('the genotype card says the call', () => {
  it('the calls of both ALT alleles’ columns are the file’s text, as one call whatever the allele', async () => {
    const t = await matrix();
    expect(t.callTexts![t.calls![0]]).toBe('0/0');
    expect(Calls.map((_, s) => t.callTexts![t.calls![s]])).toEqual(Calls);
    // the second ALT's column holds the same calls (a call is the sample's, not the allele's)
    expect(Calls.map((_, s) => t.callTexts![t.calls![t.states + s]])).toEqual(Calls);
  });

  const cases: [string, number, string][] = [
    ['homozygous reference', 0, '0/0 (homozygous reference)'],
    ['heterozygous', 1, '0/1 (heterozygous)'],
    ['homozygous alternate', 2, '1/1 (homozygous alternate)'],
    ['phased', 3, '0|1 (heterozygous)'],
    ['multi-allelic on the first ALT’s column', 4, '1/2 (heterozygous)'],
    ['haploid reference', 5, '0 (reference)'],
    ['haploid alternate', 6, '1 (alternate)'],
    ['missing', 7, './. (no call)'],
    ['partly called', 8, './1 (partly called)'],
  ];
  for (const [name, s, said] of cases)
    it(`${name}: the card reads ${said}`, async () => {
      const t = await matrix();
      const card = cardOf(t, s)!;
      expect(card.rows[1]).toEqual(['Sample', `S${s} · ${said}`]);
      expect(card.text).toContain(`Sample S${s} · ${said}`);
    });

  it('multi-allelic on each ALT’s column: the same call and word on the second ALT’s as on the first', async () => {
    const t = await matrix();
    expect(t.start.length).toBe(2);
    // both items share a position; the shortest-item rule picks the first: read the second's cell by the table
    const second = { ...t, start: Int32Array.of(t.start[1]), end: Int32Array.of(t.end[1]) };
    const cells = {
      ...second,
      codes: t.codes!.slice(t.states),
      calls: t.calls!.slice(t.states),
      id: undefined,
    } as MatrixItems;
    const card = cardOf(cells, 4)!;
    expect(card.rows[1][1]).toBe('S4 · 1/2 (heterozygous)');
  });

  it('callWords: what a call is, by its alleles', () => {
    expect(callWords('0/0')).toBe('0/0 (homozygous reference)');
    expect(callWords('2|2')).toBe('2|2 (homozygous alternate)');
    expect(callWords('0|1|1|0')).toBe('0|1|1|0 (heterozygous)');
    expect(callWords('.')).toBe('. (no call)');
  });
});
