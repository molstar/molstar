/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { columnValues } from '@molstar/core/repr/tree';
import { readTabixIndex, regionBytes } from '@molstar/io/common/tabix';
import { describe, expect, it, vi } from 'vitest';
import { sourceSuite } from '../../data/_test/source-suite.js';
import {
  DefaultOpenLimits,
  type MatrixItems,
  type MatrixSummary,
  type OpenContext,
  type RangesList,
  type RegionAnswer,
  type RegionRequest,
  type SourceEntry,
} from '../../data/protocol.js';
import { checkAnswer } from '../../data/tables.js';
import { densityOf, LevelCandidates, tabixLevels } from '../internal/levels.js';
import { genotypeClass } from '../vcf-genotypes.js';
import { alleleKey, VcfRowsBuilder } from '../vcf-parse.js';
import { VcfRuntime, vcfRuntime } from '../vcf.js';
import { memoryAssets, plainRows } from './bed-fixtures.js';
import { denseVcf, type GeneratedVcf, generateVcf, indexedVcf } from './vcf-fixtures.js';

// VCF by region past level 1 (_spec/README.md "Data", VCF: levels, the matrix): levels listed from the index's bytes (a
// dense and a sparse generated file), each tile's read cost the index's; a region read equal to the whole read clipped
// at every level (items tiles and densities counted exactly); the matrix's genotype classes and values against a
// reference reading of the same text (EC-220: uncalled 255; EC-222: a long deletion one item across the SNVs it
// covers), narrowed by `samples`, summarized above level 1, bounded before it allocates.

const signal = new AbortController().signal;
const names = ['chr1', 'chr2'];

const context = (files: Record<string, Uint8Array | string>): OpenContext => ({
  assets: memoryAssets(files),
  signal,
  limits: DefaultOpenLimits,
});
const entry = (more: Partial<SourceEntry> = {}): SourceEntry => ({
  address: '/data/v',
  format: 'vcf',
  uri: 'f',
  ...more,
});

/** A level's tile `k` on a system. */
const tile = (system: string, level: number, k: number, shapes: RegionRequest['shapes'] = ['ranges']) => ({
  system,
  start: k * level * 1024,
  end: (k + 1) * level * 1024,
  level,
  shapes,
});

/** The levels io's own `regionBytes` allows, every tile of a level checked, finest first to the first that fails. */
async function levelsOf(index: Uint8Array, extents: Record<string, number>, maxRegionBytes: number) {
  const ix = await readTabixIndex(index);
  const out: number[] = [];
  for (const level of LevelCandidates) {
    const fits = Object.entries(extents).every(([name, extent]) => {
      for (let s = 0; s < extent; s += level * 1024)
        if (level > 1 && regionBytes(ix, name, s, Math.min(extent, s + level * 1024)) > maxRegionBytes) return false;
      return true;
    });
    if (!fits) break;
    out.push(level);
  }
  return out;
}

describe('levels from the index’s bytes', () => {
  it('a sparse file lists every candidate; a dense one only those whose every tile fits `maxRegionBytes`', async () => {
    const sparse = generateVcf(400, names, 3, 6_000_000, { contigs: names });
    const s = indexedVcf(sparse);
    const a = await VcfRuntime.open(entry({ index: 'i' }), context({ f: s.file, i: s.index }));
    expect(a.coverage.levels).toEqual([1, 16, 256, 4096]);
    // dense: an SNV every 8 positions with 40 samples' calls on 300 kb, under a bound of 512 kB a tile: 16 kb tiles fit,
    // 256 kb ones do not
    const d = denseVcf(['chr1'], 300_000, 8, 40);
    const bound = 512 * 1024;
    const r = vcfRuntime({ tabix: { maxRegionBytes: bound } });
    const b = await r.open(entry({ index: 'i' }), context({ f: d.file, i: d.index }));
    const extents = { chr1: 300_000 };
    expect(b.coverage.levels).toEqual(await levelsOf(d.index, extents, bound));
    expect(b.coverage.levels).toEqual([1, 16]);
    // `maxLevel` caps them
    const capped = await vcfRuntime({ maxLevel: 256 }).open(entry({ index: 'i' }), context({ f: s.file, i: s.index }));
    expect(capped.coverage.levels).toEqual([1, 16, 256]);
  });

  it('the first level that does not fit ends the list: the coarser ones are never asked their cost', () => {
    const asked = new Set<number>();
    // a tile's cost its width: 16 kb tiles fit a bound of 20,000, 256 kb ones do not
    const cost = (_: string, a: number, b: number) => {
      asked.add(b - a);
      return b - a;
    };
    const systems = Array.from({ length: 1000 }, (_, k) => ({ name: `scaffold${k}`, extent: 8 * 1024 * 1024 }));
    expect(tabixLevels(cost, systems, { maxLevel: 4096, tileBins: 1024, maxRegionBytes: 20_000 })).toEqual([1, 16]);
    expect(asked.has(4096 * 1024)).toBe(false);
  });

  it('a sequence without `##contig` is read to the end of its last indexed 4 Mb tile', async () => {
    const g = generateVcf(300, names, 5, 9_000_000);
    const b = indexedVcf(g);
    const s = await VcfRuntime.open(entry({ index: 'i' }), context({ f: b.file, i: b.index }));
    const ix = await readTabixIndex(b.index);
    const extents = Object.fromEntries(names.map((n) => [n, 12 * 1024 * 1024]));
    expect(s.coverage.levels).toEqual(await levelsOf(b.index, extents, 8 * 1024 * 1024));
    // a tile's read cost is the index's, nothing read; past the coverage's levels none is known
    for (const level of s.coverage.levels)
      for (let k = 0; k < 4; k++) {
        const t = tile('chr1', level, k);
        expect(s.readCost!(t)).toBe(regionBytes(ix, 'chr1', t.start, t.end));
      }
    expect(s.readCost!({ ...tile('chr1', 1, 0), level: 3 })).toBeUndefined();
  });

  it('passes the source suite at every listed level, ranges and the matrix', async () => {
    const g = generateVcf(500, names, 6, 2_000_000, { contigs: names, samples: 12, genotypes: true });
    const b = indexedVcf(g);
    await sourceSuite(VcfRuntime, entry({ index: 'i' }), context({ f: b.file, i: b.index }));
    await sourceSuite(VcfRuntime, entry(), context({ f: g.text }));
  });
});

/** A whole and a by-region source over one generated file. */
async function both(g: GeneratedVcf, more: Partial<SourceEntry> = {}) {
  const b = indexedVcf(g, 2500);
  const whole = await VcfRuntime.open(entry(more), context({ f: g.text }));
  const region = await VcfRuntime.open(entry({ index: 'i', ...more }), context({ f: b.file, i: b.index }));
  return { whole, region };
}

describe('a region read equals the whole read clipped, at every level', () => {
  it('items tiles equal the whole read of the tile; densities are counted from it exactly, never estimated (random tiles, EC-171)', async () => {
    const rng = new PCG(9);
    for (let seed = 1; seed <= 3; seed++) {
      const g = generateVcf(500 + 100 * seed, names, seed, 3_000_000, { contigs: names });
      const { whole, region } = await both(g);
      for (const level of region.coverage.levels) {
        for (let n = 0; n < 6; n++) {
          const system = names[Math.floor(rng.float() * names.length)];
          const k = Math.floor((rng.float() * 3_000_000) / (level * 1024));
          const t = tile(system, level, k);
          const w = await whole.read({ ...t, level: 1 }, signal);
          const items = await region.read({ ...t, items: true }, signal);
          expect(checkAnswer(items, { ...t, items: true })).toBeUndefined();
          expect(plainRows(items), `${level} ${system} ${k}`).toEqual(plainRows(w));
          const d = await region.read(t, signal);
          expect(checkAnswer(d, t)).toBeUndefined();
          if (!d.ok || !w.ok) throw new Error('not read');
          if (level === 1) expect(plainRows(d)).toEqual(plainRows(w));
          else expect(d.tables.ranges).toEqual(densityOf(w.tables.ranges as RangesList, t));
        }
      }
      whole.dispose();
      region.dispose();
    }
  });

  it('a density counts each allele once in every bin its extent touches: a long deletion across its bins (EC-222)', async () => {
    const text = [
      '##fileformat=VCFv4.2',
      '##contig=<ID=chr1,length=2000000>',
      '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
      'chr1\t100\t.\tA\tC,G\t.\t.\t.',
      'chr1\t20000\t.\tA\t<DEL>\t.\t.\tEND=60000',
      'chr1\t30000\t.\tA\tT\t.\t.\t.',
    ].join('\n');
    const g = { text: `${text}\n` } as GeneratedVcf;
    const s = await VcfRuntime.open(entry(), context({ f: g.text }));
    const w = await s.read({ system: 'chr1', start: 0, end: 16_384 * 4, level: 1, shapes: ['ranges'] }, signal);
    if (!w.ok) throw new Error('not read');
    const d = densityOf(w.tables.ranges as RangesList, { start: 0, end: 16_384 * 4, level: 16_384 }, true);
    // bin 0: the two alleles at 100; bin 1: the deletion (from 19,999) and the SNV at 30,000; bins 2, 3: the deletion
    expect([Array.from(d.bin), Array.from(d.count)]).toEqual([
      [0, 1, 2, 3],
      [2, 2, 1, 1],
    ]);
  });
});

/** The reference reading of a genotype column's `GT` for allele `k` (1-based): its class and copies (NaN uncalled). */
function referenceCall(column: string | undefined, k: number) {
  const gt = column === undefined ? '' : column.split(':')[0];
  const alleles = gt === '' ? [] : gt.split(/[/|]/);
  const dots = alleles.filter((a) => a === '.').length;
  const copies = alleles.filter((a) => a !== '.' && Number(a) === k).length;
  return {
    code: genotypeClass(copies, alleles.length, dots),
    copies: !alleles.length || dots === alleles.length ? Number.NaN : copies,
  };
}

/** The reference matrix of a whole read's lines over `[start, end)`: rows by the sites' ids, the samples held. */
function referenceMatrix(
  g: GeneratedVcf,
  system: string,
  ids: readonly string[],
  held: readonly number[],
  field?: string,
) {
  const byId = new Map<string, number[]>();
  const counts = new Map<string, number>();
  for (const line of g.text.split('\n')) {
    if (!line || line.startsWith('#')) continue;
    const cols = line.split('\t');
    if (cols[0] !== system) continue;
    const alts = cols[4].split(',');
    const format = cols[8].split(':');
    for (let k = 0; k < alts.length; k++) {
      if (['.', '*', '<NON_REF>', '<*>'].includes(alts[k])) continue;
      // the id rule's key (its `@end` when END extends the record past REF)
      const info = cols[7];
      const end = /(^|;)END=(\d+)/.exec(info);
      const at = Number(cols[1]) - 1;
      const refEnd = at + cols[3].length;
      const endKey = end && Number(end[2]) > at && Number(end[2]) !== refEnd ? `@${end[2]}` : '';
      const key = `${cols[1]}:${alleleKey(cols[3])}>${alleleKey(alts[k])}${endKey}`;
      const n = (counts.get(key) ?? 0) + 1;
      counts.set(key, n);
      const id = n === 1 ? key : `${key}#${n}`;
      byId.set(
        id,
        held.map((s) => {
          const c = cols[9 + s];
          if (field === undefined) return referenceCall(c, k + 1).code;
          if (field === 'GT') return referenceCall(c, k + 1).copies;
          const parts = (c ?? '').split(':');
          const v = parts[format.indexOf(field)];
          const text = field === 'DS' ? v?.split(',')[k] : v;
          return text === undefined || text === '' || text === '.' ? Number.NaN : Number(text);
        }),
      );
    }
  }
  return ids.map((id) => byId.get(id));
}

describe('the matrix (genotypes, read only when asked)', () => {
  const g = generateVcf(400, names, 12, 600_000, { contigs: names, samples: 20, genotypes: true });

  it('codes the genotype class, `GT` the copies, `DS` and `GQ` their values: as a reference reading, whole and by region', async () => {
    const { whole, region } = await both(g);
    const shape = region.coverage.shapes.find((x) => x.shape === 'matrix')!;
    expect(shape.states).toEqual(Array.from({ length: 20 }, (_, k) => `S${k + 1}`));
    expect(shape.fields?.map((f) => f.name)).toEqual(['GT', 'DS', 'GQ']);
    const held = Array.from({ length: 20 }, (_, k) => k);
    for (const field of [undefined, 'GT', 'DS', 'GQ']) {
      const request: RegionRequest = {
        system: 'chr1',
        start: 0,
        end: 300_000,
        level: 1,
        shapes: ['ranges', 'matrix'],
        ...(field && { field }),
      };
      const w = await whole.read(request, signal);
      const r = await region.read(request, signal);
      expect(checkAnswer(r, request)).toBeUndefined();
      if (!w.ok || !r.ok) throw new Error('not read');
      const m = r.tables.matrix as MatrixItems;
      expect(m).toEqual(w.tables.matrix);
      const ids = columnValues(m.id!) as string[];
      expect(ids).toEqual(columnValues((r.tables.ranges as RangesList).id!));
      const got = ids.map((_, i) => Array.from((field ? m.value! : m.codes!).subarray(i * 20, (i + 1) * 20)));
      const want = referenceMatrix(g, 'chr1', ids, held, field);
      for (let i = 0; i < ids.length; i++) {
        expect(got[i].map((v) => (Number.isNaN(v) ? null : field === 'DS' ? Math.round(v * 100) / 100 : v))).toEqual(
          want[i]!.map((v) => (Number.isNaN(v) ? null : v)),
        );
      }
    }
    // the classes the generator writes, each met: none, some, every copy (a haploid's), partly called, uncalled (EC-220)
    const a = await region.read({ system: 'chr1', start: 0, end: 600_000, level: 1, shapes: ['matrix'] }, signal);
    if (!a.ok) throw new Error('not read');
    expect(new Set((a.tables.matrix as MatrixItems).codes)).toEqual(new Set([0, 1, 2, 3, 255]));
  });

  it('`samples` narrows the states, in the header’s order; names the file lacks are said once', async () => {
    const { region } = await both(g, { samples: ['S7', 'S2', 'nobody'] } as Partial<SourceEntry>);
    expect(region.coverage.shapes.find((x) => x.shape === 'matrix')?.states).toEqual(['S2', 'S7']);
    expect(region.issues?.find((i) => i.code === 'seq.source.vcf.samples')?.message).toMatch(/nobody/);
    const request: RegionRequest = { system: 'chr1', start: 0, end: 600_000, level: 1, shapes: ['matrix'] };
    const a = await region.read(request, signal);
    if (!a.ok) throw new Error('not read');
    const m = a.tables.matrix as MatrixItems;
    const ids = columnValues(m.id!) as string[];
    const want = referenceMatrix(g, 'chr1', ids, [1, 6]);
    for (let i = 0; i < ids.length; i++) expect(Array.from(m.codes!.subarray(2 * i, 2 * i + 2))).toEqual(want[i]);
  });

  it('a file without samples serves no matrix; a matrix past `maxMatrixBytes` is too large, refused before it allocates', async () => {
    const sites = generateVcf(100, names, 2, 100_000);
    const s = await VcfRuntime.open(entry(), context({ f: sites.text }));
    expect(s.coverage.shapes.some((x) => x.shape === 'matrix')).toBe(false);
    const r = await s.read({ system: 'chr1', start: 0, end: 1000, level: 1, shapes: ['matrix'] }, signal);
    expect(r).toMatchObject({ ok: false, kind: 'not-served' });
    const small = await vcfRuntime({ maxMatrixBytes: 100 }).open(entry(), context({ f: g.text }));
    const big = await small.read({ system: 'chr1', start: 0, end: 600_000, level: 1, shapes: ['matrix'] }, signal);
    expect(big).toMatchObject({ ok: false, kind: 'too-large' });
  });

  it('a matrix read’s line buffers count in `maxMatrixBytes` (states × the highest allele), never held past the read', async () => {
    const r: RegionRequest = { system: 'chr1', start: 0, end: 600_000, level: 1, shapes: ['matrix'] };
    const all = (await (await VcfRuntime.open(entry(), context({ f: g.text }))).read(r, signal)) as RegionAnswer & {
      ok: true;
    };
    const cells = (all.tables.matrix as MatrixItems).codes!.length;
    // the items' cells (a code and a call each) alone fit, with the buffers they do not
    const s = await vcfRuntime({ maxMatrixBytes: 2 * cells }).open(entry(), context({ f: g.text }));
    expect(await s.read(r, signal)).toMatchObject({ ok: false, kind: 'too-large' });
    const t = await vcfRuntime({ maxMatrixBytes: 2 * cells + 20 * (3 + 10) }).open(entry(), context({ f: g.text }));
    expect((await t.read(r, signal)).ok).toBe(true);
  });

  it('a DS or GQ value read from its bytes as `Number` reads its text: decimals, a sign, an exponent, `.` missing', async () => {
    const text = [
      '##fileformat=VCFv4.2',
      '##contig=<ID=chr1,length=1000>',
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
      '##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">',
      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\tE\tF\tG',
      'chr1\t10\t.\tA\tC\t.\t.\t.\tGT:DS\t0/1:0.25\t0/1:-1.5e-3\t0/1:.\t0/1:12\t0/1:+.5\t0/1:Inf\t0/1:1x',
    ].join('\n');
    const s = await VcfRuntime.open(entry(), context({ f: `${text}\n` }));
    const a = (await s.read(
      { system: 'chr1', start: 0, end: 100, level: 1, shapes: ['matrix'], field: 'DS' },
      signal,
    )) as RegionAnswer & { ok: true };
    const want = [0.25, -1.5e-3, Number.NaN, 12, 0.5, Number.POSITIVE_INFINITY, Number.NaN].map(Math.fround);
    expect(Array.from((a.tables.matrix as MatrixItems).value!)).toEqual(want);
  });

  it('a whole read keeps the file’s bytes only when it holds samples (a sites file’s rows hold no copy of it)', async () => {
    const kept: number[] = [];
    const finish = VcfRowsBuilder.prototype.finish;
    const spy = vi.spyOn(VcfRowsBuilder.prototype, 'finish').mockImplementation(function (this: VcfRowsBuilder, b) {
      kept.push(b.length);
      return finish.call(this, b);
    });
    await VcfRuntime.open(entry(), context({ f: generateVcf(50, names, 2, 100_000).text }));
    expect(kept.length).toBeGreaterThan(0);
    expect(kept.every((n) => n === 0)).toBe(true);
    kept.length = 0;
    await VcfRuntime.open(entry(), context({ f: g.text }));
    expect(kept.every((n) => n > 0)).toBe(true);
    spy.mockRestore();
  });

  it('a long deletion is one item across the SNVs it covers, each sample’s class its own (EC-222)', async () => {
    const text = [
      '##fileformat=VCFv4.2',
      '##contig=<ID=chr1,length=100000>',
      '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD',
      'chr1\t1000\t.\tA\t<DEL>\t.\t.\tEND=5000\tGT\t0/1\t0/0\t1\t./.',
      'chr1\t2000\t.\tC\tT\t.\t.\t.\tGT\t0/0\t1/1\t0\t./1',
      'chr1\t3000\t.\tG\tA,T\t.\t.\t.\tGT\t1/2\t2/2\t.\t0|1',
    ].join('\n');
    const s = await VcfRuntime.open(entry(), context({ f: `${text}\n` }));
    const a = await s.read({ system: 'chr1', start: 2500, end: 2600, level: 1, shapes: ['matrix'] }, signal);
    if (!a.ok) throw new Error('not read');
    const m = a.tables.matrix as MatrixItems;
    // the region inside the deletion: the deletion alone, from its start to its END
    expect([Array.from(m.start), Array.from(m.end), Array.from(m.codes!)]).toEqual([[999], [5000], [1, 0, 2, 255]]);
    const all = await s.read(
      { system: 'chr1', start: 0, end: 10_000, level: 1, shapes: ['matrix'], field: 'GT' },
      signal,
    );
    if (!all.ok) throw new Error('not read');
    const v = all.tables.matrix as MatrixItems;
    expect(columnValues(v.id!)).toEqual(['1000:A><DEL>@5000', '2000:C>T', '3000:G>A', '3000:G>T']);
    expect(Array.from(v.value!, (x) => (Number.isNaN(x) ? null : x))).toEqual([
      ...[1, 0, 1, null],
      ...[0, 2, 0, 1],
      ...[1, 0, null, 1],
      ...[1, 2, null, 0],
    ]);
  });

  it('above level 1, a summary: each bin’s majority class and share, or mean, min, max and count, counted from the tile', async () => {
    const { whole, region } = await both(g);
    const level = region.coverage.levels[1];
    const t = { ...tile('chr1', level, 0, ['matrix']) };
    const a = await region.read(t, signal);
    expect(checkAnswer(a, t)).toBeUndefined();
    const w = (await whole.read({ ...t, level: 1 }, signal)) as RegionAnswer & { ok: true };
    const items = w.tables.matrix as MatrixItems;
    const s = (a as RegionAnswer & { ok: true }).tables.matrix as MatrixSummary;
    expect(s.form).toBe('summary');
    // bin by bin, state by state: the reference majority of the items touching the bin
    const S = 20;
    for (let bin = 0; bin < 1024; bin += 97)
      for (let st = 0; st < S; st++) {
        const counts = [0, 0, 0, 0];
        for (let k = 0; k < items.start.length; k++) {
          const b0 = Math.floor(Math.max(items.start[k], t.start) / level);
          const b1 = Math.floor((Math.max(items.end[k], items.start[k] + 1) - 1) / level);
          if (bin < b0 || bin > b1) continue;
          const c = items.codes![k * S + st];
          if (c <= 3) counts[c]++;
        }
        const called = counts.reduce((x, y) => x + y, 0);
        const best = counts.indexOf(Math.max(...counts));
        expect([s.classes![bin * S + st], s.shares![bin * S + st]]).toEqual(
          called ? [best, Math.round((255 * counts[best]) / called)] : [255, 0],
        );
      }
    const gq = await region.read({ ...t, field: 'GQ' }, signal);
    expect(checkAnswer(gq, { ...t, field: 'GQ' })).toBeUndefined();
    expect(((gq as RegionAnswer & { ok: true }).tables.matrix as MatrixSummary).mean).toBeInstanceOf(Float32Array);
    // a matrix asked as items at a coarse level: the items of that tile, as at level 1
    const it = await region.read({ ...t, items: true }, signal);
    expect((it as RegionAnswer & { ok: true }).tables.matrix).toEqual(w.tables.matrix);
  });

  it('fuzzed genotype columns read without a throw: every byte of a sample column changed, a column cut short (EC-40)', async () => {
    const rng = new PCG(31);
    const text = new TextEncoder().encode(generateVcf(60, names, 8, 100_000, { samples: 8, genotypes: true }).text);
    const body = text.indexOf(10, text.lastIndexOf(35)) + 1;
    for (let k = 0; k < 40; k++) {
      const bytes = text.slice(0, body + Math.floor(rng.float() * (text.length - body)));
      for (let n = 0; n < 50; n++) {
        const at = body + Math.floor(rng.float() * (bytes.length - body));
        bytes[at] = ['/', '|', '.', ':', ',', '\t', '9', 'x'][Math.floor(rng.float() * 8)].charCodeAt(0);
      }
      const s = await VcfRuntime.open(entry(), context({ f: bytes }));
      for (const field of [undefined, 'GT', 'DS', 'GQ']) {
        const r: RegionRequest = {
          system: 'chr1',
          start: 0,
          end: 100_000,
          level: 1,
          shapes: ['matrix'],
          ...(field && { field }),
        };
        const a = await s.read(r, signal);
        expect(checkAnswer(a, r)).toBeUndefined();
      }
    }
  });
});
