/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { gzipSync } from 'node:zlib';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { columnValues } from '@molstar/core/repr/tree';
import { bgzfReader } from '@molstar/io/common/bgzf';
import { readRegion, readTabixIndex } from '@molstar/io/common/tabix';
import { bgzip, memoryRanges } from '@molstar/test-harness/genome/bgzf-writer';
import { describe, expect, it } from 'vitest';
import { sourceSuite } from '../../data/_test/source-suite.js';
import { delayedRuntime } from '../../data/async-source.js';
import {
  DefaultOpenLimits,
  type OpenContext,
  type OpenedSource,
  type RegionAnswer,
  type SourceEntry,
} from '../../data/protocol.js';
import { SourceOpenError } from '../../data/seq-runtime.js';
import { checkAnswer } from '../../data/tables.js';
import { DefaultVcfOptions, MissingInt, readVcfHeader, VcfRuntime, vcfRuntime } from '../vcf.js';
import { alleleCategory, alleleKey, vcfExtent } from '../vcf-parse.js';
import { memoryAssets, plainRows } from './bed-fixtures.js';
import { generateVcf, indexedVcf, VcfHeaderLines } from './vcf-fixtures.js';

// VCF sites (_spec/README.md "Data", VCF): the source suite on generated files, whole and by region; a region read
// equal to the whole read clipped (a property over generated files and random regions, and over tilings, rows deduped
// by id); each line's extent equal to io's (htslib's `tbx_parse1`); ids, alleles, INFO's types and `Number`s; fuzzed and
// truncated files without a throw; every bound refused before it allocates. The recorded ClinVar slice is in
// `vcf-clinvar.test.ts`.

const signal = new AbortController().signal;
const names = ['chr1', 'chr2', 'chrX'];

function context(files: Record<string, Uint8Array | string>, ranges = true, limits = DefaultOpenLimits): OpenContext {
  return { assets: memoryAssets(files, ranges), signal, limits };
}

const entry = (more: Partial<SourceEntry> = {}): SourceEntry => ({
  address: '/data/v',
  format: 'vcf',
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

/** An answer's rows as records of their columns (`plainRows`), or why there are none. */
function rowsOf(a: RegionAnswer): Record<string, unknown>[] {
  return plainRows(a) as Record<string, unknown>[];
}

const vcf = (...lines: string[]) =>
  `${[...VcfHeaderLines, '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', ...lines].join('\n')}\n`;

async function readAll(text: string, system = 'chr1', more: Partial<SourceEntry> = {}) {
  const s = await VcfRuntime.open(entry(more), context({ f: text }));
  const a = await s.read({ system, start: 0, end: 1_000_000, level: 1, shapes: ['ranges'] }, signal);
  return { s, rows: rowsOf(a) };
}

describe('the source suite', () => {
  const g = generateVcf(600, names, 4, 400_000, { contigs: ['chr1', 'chr2'] });

  it('read whole: plain, gzip and BGZF by their magic, with and without range reads', async () => {
    const text = new TextEncoder().encode(g.text);
    for (const bytes of [text, new Uint8Array(gzipSync(text)), bgzip(text, 3000).bytes])
      for (const ranges of [true, false]) await sourceSuite(VcfRuntime, entry(), context({ f: bytes }, ranges));
  });

  it('by region through tabix (level 1), and the synthetic twin delayed', async () => {
    const b = indexedVcf(g);
    await sourceSuite(VcfRuntime, entry({ index: 'i' }), context({ f: b.file, i: b.index }));
    await sourceSuite(VcfRuntime, entry({ index: 'i' }), context({ f: b.file, i: b.index }, false));
    await sourceSuite(delayedRuntime(VcfRuntime, { delayMs: 1 }), entry(), context({ f: g.text }));
  });
});

describe('a region read equals the whole read clipped', () => {
  it('generated files, random regions, every answer equal row for row', async () => {
    const rng = new PCG(11);
    for (let seed = 1; seed <= 6; seed++) {
      const g = generateVcf(300 + seed * 53, names, seed, 400_000, { contigs: seed % 2 ? ['chr1'] : [] });
      const { file, index } = indexedVcf(g, 1000 + seed * 500);
      const whole = await VcfRuntime.open(entry(), context({ f: g.text }));
      const region = await VcfRuntime.open(entry({ index: 'i' }), context({ f: file, i: index }));
      expect(region.coverage).toEqual({ ...whole.coverage, levels: [1] });
      for (let k = 0; k < 25; k++) {
        const at = g.lines[Math.floor(rng.float() * g.lines.length)];
        // every fifth region starts at a line's own start, every fifth a position long
        const start = k % 5 === 1 ? at.begin : Math.floor(rng.float() * 400_000);
        const end = start + Math.floor(rng.float() * (k % 5 === 0 ? 2 : 30_000));
        const system = k % 5 === 1 ? names[at.seq] : names[Math.floor(rng.float() * names.length)];
        const request = { system, start, end, level: 1, shapes: ['ranges' as const] };
        const w = await whole.read(request, signal);
        const r = await region.read(request, signal);
        expect(checkAnswer(r, request)).toBeUndefined();
        expect(rowsOf(r), `${seed} ${system}:${start}-${end}`).toEqual(rowsOf(w));
      }
      whole.dispose();
      region.dispose();
    }
  });

  it('a tiling of a sequence, its rows deduped by id, is the whole read (a long deletion one item across its tiles, EC-222)', async () => {
    const g = generateVcf(500, names, 7, 200_000);
    const { file, index } = indexedVcf(g, 1500);
    const whole = await VcfRuntime.open(entry(), context({ f: g.text }));
    const region = await VcfRuntime.open(entry({ index: 'i' }), context({ f: file, i: index }));
    for (const tile of [1024, 8192, 65_536]) {
      const byId = new Map<unknown, Record<string, unknown>>();
      let longest = 0;
      for (let start = 0; start < 220_000; start += tile) {
        const a = await region.read({ system: 'chr2', start, end: start + tile, level: 1, shapes: ['ranges'] }, signal);
        for (const row of rowsOf(a)) {
          const seen = byId.get(row.id);
          // a row meeting several tiles is the same row in each: one item, its whole extent
          if (seen) expect(row).toEqual(seen);
          else byId.set(row.id, row);
          longest = Math.max(longest, (row.end as number) - (row.start as number));
        }
      }
      const all = rowsOf(
        await whole.read({ system: 'chr2', start: 0, end: 220_000, level: 1, shapes: ['ranges'] }, signal),
      );
      expect([...byId.values()].sort((a, b) => (a.start as number) - (b.start as number))).toEqual(
        [...all].sort((a, b) => (a.start as number) - (b.start as number)),
      );
      if (tile < 20_000) expect(longest).toBeGreaterThan(tile);
    }
  });
});

describe('a line’s extent, by the index’s rule (htslib’s tbx_parse1)', () => {
  it('equals io’s for every generated line: `END=` on an SNV line, a symbolic deletion, a breakend', async () => {
    const g = generateVcf(800, names, 3);
    const { file, index } = indexedVcf(g, 2000);
    const ix = await readTabixIndex(index);
    const reader = bgzfReader(memoryRanges(file));
    let checked = 0;
    for (const name of names) {
      const l = await readRegion(reader, ix, name, 0, 2 ** 29);
      for (let i = 0; i < l.lineStart.length; i++) {
        const cells = new TextDecoder().decode(l.bytes.subarray(l.lineStart[i], l.lineEnd[i])).split('\t');
        expect(vcfExtent(Number(cells[1]), cells[3], cells[7]), cells.join(' ')).toEqual([l.begin[i], l.end[i]]);
        checked++;
      }
    }
    expect(checked).toBe(g.lines.length);
    // the named cases
    expect(vcfExtent(100, 'A', 'DP=3;END=250')).toEqual([99, 250]);
    expect(vcfExtent(100, 'A', 'END=250;DP=3')).toEqual([99, 250]);
    expect(vcfExtent(100, 'N', 'SVTYPE=DEL;END=5100;SVLEN=-5000')).toEqual([99, 5100]);
    expect(vcfExtent(100, 'G', 'SVTYPE=BND;MATEID=b2')).toEqual([99, 100]);
    expect(vcfExtent(100, 'ACGT', 'SVEND=900')).toEqual([99, 103]);
    expect(vcfExtent(100, 'A', 'END=99')).toEqual([99, 100]);
    expect(vcfExtent(100, 'A', 'END=.')).toEqual([99, 100]);
  });
});

describe('rows: a row an ALT allele, its id, its fields (EC-196)', () => {
  it('a multi-allelic line is a row an allele, each with the line’s extent; ALT `.`, `*` and `<NON_REF>` no item', async () => {
    const { s, rows } = await readAll(
      vcf(
        'chr1\t100\trs1\tA\tC,G,T\t50\tPASS\tAF=0.1,0.2,0.3;RC=10,1,2,3',
        'chr1\t200\t.\tC\tT,*\t.\t.\t.',
        'chr1\t300\t.\tG\t.\t.\t.\t.',
        'chr1\t400\t.\tG\t<NON_REF>\t.\t.\t.',
      ),
    );
    expect(rows.map((r) => [r.id, r.start, r.end, r.alt, r.AF, r.RC])).toEqual([
      ['100:A>C', 99, 100, 'C', expect.closeTo(0.1, 6), '10,1'],
      ['100:A>G', 99, 100, 'G', expect.closeTo(0.2, 6), '10,2'],
      ['100:A>T', 99, 100, 'T', expect.closeTo(0.3, 6), '10,3'],
      ['200:C>T', 199, 200, 'T', null, null],
    ]);
    expect(s.issues?.map((i) => i.code)).toContain('seq.source.vcf.allele');
    expect(s.issues?.find((i) => i.code === 'seq.source.vcf.allele')?.message).toMatch(/line 1[34]/);
  });

  it('ids: `pos:ref>alt`, a repeated key `#2` in file order, alleles past 32 characters `<length>:<FNV-1a hex>`', async () => {
    const long = 'ACGT'.repeat(10);
    const { rows } = await readAll(
      vcf(
        'chr1\t100\trs1\tA\tC\t.\t.\t.',
        'chr1\t100\trs1b\tA\tC\t.\t.\t.',
        'chr1\t100\t.\tA\tC\t.\t.\t.',
        `chr1\t500\tsv1\t${long}\tA\t.\t.\t.`,
        `chr1\t900\tins1\tA\tA${long}\t.\t.\t.`,
      ),
    );
    expect(rows.map((r) => r.id)).toEqual([
      '100:A>C',
      '100:A>C#2',
      '100:A>C#3',
      `500:${alleleKey(long)}>A`,
      `900:A>${alleleKey(`A${long}`)}`,
    ]);
    expect(alleleKey(long)).toMatch(/^40:[0-9a-f]{8}$/);
    expect(alleleKey(`A${long}`)).toMatch(/^41:[0-9a-f]{8}$/);
    expect(alleleKey(long)).not.toBe(alleleKey(`${long.slice(1)}A`));
    // the VCF ID the name and the label, the ref and alt whole
    expect(rows.map((r) => [r.name, r.label])).toEqual([
      ['rs1', 'rs1'],
      ['rs1b', 'rs1b'],
      [null, null],
      ['sv1', 'sv1'],
      ['ins1', 'ins1'],
    ]);
    expect(rows[3].ref).toBe(long);
    expect(rows[3].end).toBe(499 + 40);
  });

  it('categories: SNV, MNV, insertion, deletion, indel, structural, breakend', () => {
    expect(
      [
        ['A', 'C'],
        ['AC', 'GT'],
        ['A', 'ACGT'],
        ['ACGT', 'A'],
        ['AC', 'GTT'],
        ['A', '<DEL>'],
        ['A', '<INS:ME:ALU>'],
        ['G', 'G]17:198982]'],
        ['G', ']13:123456]T'],
        ['G', '.G'],
        ['G', 'G.'],
      ].map(([r, a]) => alleleCategory(r, a)),
    ).toEqual([
      'SNV',
      'MNV',
      'insertion',
      'deletion',
      'indel',
      'structural',
      'structural',
      'breakend',
      'breakend',
      'breakend',
      'breakend',
    ]);
  });

  it('a long deletion is one item across the SNVs it covers (EC-222, sites): its extent to its END in every region it meets', async () => {
    const text = vcf(
      'chr1\t1000\tdel1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=6000',
      'chr1\t2000\ts1\tA\tC\t.\tPASS\t.',
      'chr1\t3000\ts2\tG\tT\t.\tPASS\t.',
    );
    const s = await VcfRuntime.open(entry(), context({ f: text }));
    for (const [start, end] of [
      [1500, 2500],
      [2900, 3100],
      [5000, 5001],
    ]) {
      const rows = rowsOf(await s.read({ system: 'chr1', start, end, level: 1, shapes: ['ranges'] }, signal));
      expect(rows[0]).toMatchObject({ id: '1000:N><DEL>', start: 999, end: 6000, category: 'structural' });
    }
  });

  it('fixed fields: `qual` float32 (`.` NaN), `filter` a category (`.` none)', async () => {
    const { s, rows } = await readAll(vcf('chr1\t10\t.\tA\tC\t29.5\tq10;s50\t.', 'chr1\t20\t.\tA\tC\t.\t.\t.'));
    expect(rows.map((r) => [r.qual, r.filter])).toEqual([
      [29.5, 'q10;s50'],
      [null, null],
    ]);
    expect(s.coverage.shapes[0].fields?.slice(0, 7)).toEqual([
      { name: 'label', type: 'category' },
      { name: 'name', type: 'category' },
      { name: 'ref', type: 'category' },
      { name: 'alt', type: 'category' },
      { name: 'qual', type: 'float32' },
      { name: 'filter', type: 'category' },
      { name: 'category', type: 'category' },
    ]);
  });
});

describe('INFO: typed by the header, each `Number`', () => {
  it('`Integer` int32 (missing BCF’s), `Float` float32, `Flag` and `String` categories; `A` by allele, `R` the reference’s and the allele’s, others their text', async () => {
    const { s, rows } = await readAll(
      vcf(
        'chr1\t100\t.\tA\tC,G\t.\t.\tDP=12;AF=0.25,0.5;RC=7,3,1;DB;CLNSIG=Pathogenic,Benign;GL=-1,-2,-3',
        'chr1\t200\t.\tA\tC\t.\t.\tDP=.;AF=x',
      ),
    );
    const fields = Object.fromEntries((s.coverage.shapes[0].fields ?? []).map((f) => [f.name, f.type]));
    expect(fields).toMatchObject({
      DP: 'int32',
      AF: 'float32',
      RC: 'category',
      DB: 'category',
      CLNSIG: 'category',
      SVTYPE: 'category',
      END: 'int32',
      GL: 'category',
    });
    expect(rows.map((r) => [r.DP, r.AF, r.RC, r.DB, r.CLNSIG, r.GL])).toEqual([
      [12, 0.25, '7,3', 'true', 'Pathogenic,Benign', '-1,-2,-3'],
      [12, 0.5, '7,1', 'true', 'Pathogenic,Benign', '-1,-2,-3'],
      [MissingInt, null, null, null, null, null],
    ]);
    expect(s.issues?.find((i) => i.code === 'seq.source.vcf.cell')?.message).toMatch(/no number.*line 13/);
  });

  it('`fields` and `types` first and by their names (EC-190); the header’s order; a key undeclared said; one named as a fixed field not read', async () => {
    const text = vcf('chr1\t100\t.\tA\tC\t.\t.\tDP=12;AF=0.5;XX=1;qual=3');
    const s = await VcfRuntime.open(
      entry({ fields: { frequency: 'AF', ref: 'DP' }, types: { GL: 'string', XY: 'float32' } }),
      context({ f: text }),
    );
    const fields = s.coverage.shapes[0].fields!.slice(7).map((f) => `${f.name}:${f.type}`);
    expect(fields).toEqual([
      'frequency:float32',
      'GL:category',
      'XY:float32',
      'RC:category',
      'DB:category',
      'CLNSIG:category',
      'SVTYPE:category',
      'END:int32',
    ]);
    const codes = s.issues?.map((i) => `${i.code} ${i.message}`) ?? [];
    expect(codes.some((c) => /vcf\.cell An INFO key named as a fixed field.*DP/.test(c))).toBe(true);
    expect(codes.some((c) => /vcf\.info An INFO key the header does not declare: 2 \(XX, qual\)/.test(c))).toBe(true);
  });

  it('past `maxInfoFields` the rest are said once with their count', async () => {
    const s = await vcfRuntime({ maxInfoFields: 3 }).open(entry(), context({ f: vcf('chr1\t1\t.\tA\tC\t.\t.\t.') }));
    expect(s.coverage.shapes[0].fields?.slice(7).map((f) => f.name)).toEqual(['DP', 'AF', 'RC']);
    expect(s.issues?.find((i) => i.code === 'seq.source.vcf.info')?.message).toMatch(/past the 3 read: 5 not fields/);
  });
});

describe('the header', () => {
  it('`##fileformat`, `##contig` (names, lengths: EC-188), `##INFO`, `##FORMAT`, the samples; quoted values unescaped', () => {
    const h = readVcfHeader(
      [
        '##fileformat=VCFv4.3',
        '##contig=<ID=chr1,length=248956422,assembly=GRCh38>',
        '##contig=<ID=chrUn>',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="A, \\"b\\"",Source="x">',
        '##INFO=<ID=BAD,Number=1,Type=Wrong>',
        '##INFO=<Number=1,Type=Integer>',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2',
        'chr1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\t1/1',
      ].join('\n'),
    );
    expect(h).toEqual({
      ok: true,
      header: {
        fileformat: 'VCFv4.3',
        contigs: [{ name: 'chr1', length: 248956422 }, { name: 'chrUn' }],
        info: [{ id: 'AF', number: 'A', type: 'Float', description: 'A, "b"' }],
        format: [{ id: 'GT', number: '1', type: 'String', description: 'Genotype' }],
        samples: ['NA1', 'NA2'],
        lines: 8,
        bytes: expect.any(Number),
        unread: 2,
      },
    });
  });

  it('without `#CHROM` malformed; past `maxHeaderBytes` or `maxSamples` malformed, counted before the names are split', () => {
    expect(readVcfHeader('##fileformat=VCFv4.2\nchr1\t1\n')).toMatchObject({ ok: false, kind: 'malformed' });
    expect(readVcfHeader('##fileformat=VCFv4.2\n#CHROM\tPOS\n')).toMatchObject({ ok: false });
    const big = `##x=${'y'.repeat(2000)}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n`;
    expect(readVcfHeader(big, { maxHeaderBytes: 1000 })).toMatchObject({ ok: false, message: /over 1000 bytes/ });
    const many = `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT${'\tS'.repeat(500)}\n`;
    expect(readVcfHeader(many, { maxSamples: 499 })).toMatchObject({ ok: false, message: /500 samples/ });
    expect(readVcfHeader(many, { maxSamples: 500 })).toMatchObject({ ok: true });
  });

  it('a contig’s length is its system’s (EC-188); a row past it clipped, one starting past it left out, counted once (EC-51)', async () => {
    const text = vcf(
      'chr1\t100\t.\tA\tC\t.\t.\t.',
      'chr1\t990\t.\tACGTACGTACGTACGT\tA\t.\t.\t.',
      'chr1\t2000\t.\tA\tC\t.\t.\t.',
    ).replace(
      '##fileformat=VCFv4.2\n',
      '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n##contig=<ID=chr9,length=500>\n',
    );
    const { s, rows } = await readAll(text);
    expect(s.coverage.systems).toEqual([
      { name: 'chr1', length: 1000 },
      { name: 'chr9', length: 500 },
    ]);
    expect(rows.map((r) => [r.start, r.end])).toEqual([
      [99, 100],
      [989, 1000],
    ]);
    expect(s.issues?.find((i) => i.code === 'seq.source.vcf.range')?.message).toMatch(/: 2 \(line 15, line 16\)/);
  });

  it('a contig longer than 2^31 − 1 is refused by name, the others read', async () => {
    const text = vcf('chr1\t100\t.\tA\tC\t.\t.\t.', 'chrBig\t100\t.\tA\tC\t.\t.\t.').replace(
      '##fileformat=VCFv4.2\n',
      '##fileformat=VCFv4.2\n##contig=<ID=chrBig,length=3000000000>\n',
    );
    const { s, rows } = await readAll(text);
    expect(rows.length).toBe(1);
    expect(s.coverage.systems.map((x) => x.name)).toEqual(['chr1']);
    expect(s.issues?.find((i) => i.code === 'seq.source.too_large')?.message).toMatch(/chrBig/);
  });
});

describe('the source', () => {
  it('`base` and `ends` not read, said once (EC-195); `samples` kept and said as not yet read, never an error', async () => {
    const s = await VcfRuntime.open(
      entry({ base: 0, ends: 'exclusive', samples: ['NA1'] } as Partial<SourceEntry>),
      context({ f: vcf('chr1\t100\t.\tA\tC\t.\t.\t.') }),
    );
    expect(s.issues?.map((i) => [i.code, i.severity])).toEqual([
      ['seq.source.vcf.transform', 'info'],
      ['seq.source.vcf.samples', 'info'],
    ]);
    const a = await s.read({ system: 'chr1', start: 0, end: 200, level: 1, shapes: ['ranges'] }, signal);
    expect(rowsOf(a)[0]).toMatchObject({ start: 99, end: 100 });
  });

  it('the index the source names is read, never `<uri>.tbi` (EC-186); a missing one refuses the source', async () => {
    const g = generateVcf(50, names, 2);
    const b = indexedVcf(g);
    const asked: string[] = [];
    const inner = memoryAssets({ 'data/x.vcf.gz': b.file, 'idx/elsewhere.tbi': b.index });
    const assets = {
      ...inner,
      readBytes: (uri: string, s: AbortSignal) => {
        asked.push(uri);
        return inner.readBytes(uri, s);
      },
    };
    const s = await VcfRuntime.open(entry({ uri: 'data/x.vcf.gz', index: 'idx/elsewhere.tbi' }), {
      assets,
      signal,
      limits: DefaultOpenLimits,
    });
    expect(asked).toEqual(['idx/elsewhere.tbi']);
    s.dispose();
    expect(
      await opened(
        VcfRuntime.open(entry({ uri: 'data/x.vcf.gz', index: 'nowhere' }), {
          assets,
          signal,
          limits: DefaultOpenLimits,
        }),
      ),
    ).toEqual(['seq.source.unreadable']);
  });

  it('a line that does not parse is left out and said with its line number', async () => {
    const { s, rows } = await readAll(
      vcf(
        'chr1\t100\t.\tA\tC\t.\t.\t.',
        'chr1\tten\t.\tA\tC\t.\t.\t.',
        'chr1\t300\t.\t\tC\t.\t.\t.',
        'chr1\t400\t.\tA',
        'chr1\t500\t.\tA\tC\t.\t.\tEND=9999999999',
        'chr1\t600\t.\tA\tC\tnoqual\t.\t.',
      ),
    );
    expect(rows.map((r) => r.start)).toEqual([99, 599]);
    const line = s.issues?.filter((i) => i.code === 'seq.source.vcf.line').map((i) => i.message);
    expect(line).toEqual([
      'A position that is not a whole number to 2^31 − 1: 1 (line 13)',
      'No REF: 1 (line 14)',
      'Fewer than eight columns: 1 (line 15)',
      'An end past 2^31 − 1: 1 (line 16)',
    ]);
  });
});

describe('bounds, each refused before it allocates', () => {
  it('a record past `maxAlleles` left out (counted from its commas), said', async () => {
    const alts = Array.from({ length: 300 }, (_, i) => `A${'C'.repeat(i % 5)}T`).join(',');
    const { s, rows } = await readAll(vcf(`chr1\t100\t.\tA\t${alts}\t.\t.\t.`, 'chr1\t200\t.\tA\tC\t.\t.\t.'));
    expect(rows.length).toBe(1);
    expect(s.issues?.find((i) => i.code === 'seq.source.vcf.alleles')?.message).toMatch(/More than 256.*line 12/);
    expect(DefaultVcfOptions.maxAlleles).toBe(256);
  });

  it('a header past `maxHeaderBytes` or `maxSamples`, a line past `maxLineBytes`: malformed, whole and by region', async () => {
    const g = generateVcf(40, names, 5, 100_000, { samples: 30 });
    const b = indexedVcf(g);
    for (const [over, code] of [
      [{ maxHeaderBytes: 200 }, 'seq.source.malformed'],
      [{ maxSamples: 29 }, 'seq.source.malformed'],
      [{ maxLineBytes: 60 }, 'seq.source.malformed'],
    ] as const) {
      const r = vcfRuntime(over);
      expect(await opened(r.open(entry(), context({ f: g.text })))).toEqual([code]);
      const byRegion = await opened(r.open(entry({ index: 'i' }), context({ f: b.file, i: b.index })));
      if (byRegion[0] === 'opened') {
        // a line past the bound fails the region that reads it
        const s = await r.open(entry({ index: 'i' }), context({ f: b.file, i: b.index }));
        const a = await s.read(
          { system: names[g.lines[0].seq], start: 0, end: 100_000, level: 1, shapes: ['ranges'] },
          signal,
        );
        expect(a).toMatchObject({ ok: false, kind: 'malformed' });
      } else expect(byRegion).toEqual([code]);
    }
  });

  it('read whole past `limits.maxBytes` or `limits.maxRows`: too large, hinting an index', async () => {
    const g = generateVcf(200, names, 6);
    const small = { maxRows: 50, maxBytes: DefaultOpenLimits.maxBytes };
    expect(await opened(VcfRuntime.open(entry(), context({ f: g.text }, true, small)))).toEqual([
      'seq.source.too_large',
    ]);
    const tiny = { maxRows: DefaultOpenLimits.maxRows, maxBytes: 1000 };
    expect(await opened(VcfRuntime.open(entry(), context({ f: g.text }, true, tiny)))).toEqual([
      'seq.source.too_large',
    ]);
  });
});

describe('fuzzed and truncated files read without a throw (EC-40, EC-92)', () => {
  it('bytes flipped, lines cut, files truncated: opened or refused with an issue, every answer checked', async () => {
    const g = generateVcf(200, names, 9);
    const text = new TextEncoder().encode(g.text);
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
      variants.push(v);
    }
    for (const v of variants) {
      const got = await opened(VcfRuntime.open(entry(), context({ f: v })));
      expect(got.length).toBeGreaterThan(0);
      if (got[0] !== 'opened') continue;
      const s = await VcfRuntime.open(entry(), context({ f: v }));
      for (const sh of s.coverage.shapes.slice(0, 3)) {
        const request = { system: sh.system, start: 0, end: 500_000, level: 1, shapes: [sh.shape] };
        expect(checkAnswer(await s.read(request, signal), request)).toBeUndefined();
      }
    }
  });

  it('by region: a damaged file or index, a truncated file: refused at open or a failed region, never a throw', async () => {
    const g = generateVcf(300, names, 10);
    const { file, index } = indexedVcf(g, 1500);
    const rng = new PCG(5);
    for (let k = 0; k < 30; k++) {
      const f = k % 3 === 0 ? file.slice(0, Math.floor(rng.float() * file.length)) : file.slice();
      const i = k % 3 === 1 ? index.slice() : index;
      if (k % 3 === 1) for (let n = 0; n < 5; n++) i[Math.floor(rng.float() * i.length)] ^= 0xff;
      if (k % 3 === 2) for (let n = 0; n < 20; n++) f[Math.floor(rng.float() * f.length)] ^= 0x5a;
      const got = await opened(VcfRuntime.open(entry({ index: 'i' }), context({ f, i })));
      if (got[0] !== 'opened') continue;
      const s = await VcfRuntime.open(entry({ index: 'i' }), context({ f, i }));
      for (const sh of s.coverage.shapes) {
        const request = { system: sh.system, start: 0, end: 400_000, level: 1, shapes: [sh.shape] };
        expect(checkAnswer(await s.read(request, signal), request)).toBeUndefined();
      }
      s.dispose();
    }
  });
});

describe('performance (quick representative)', () => {
  it('a region of 1,000 records decoded (sites)', async () => {
    const g = generateVcf(1000, ['chr1'], 12, 100_000);
    const { file, index } = indexedVcf(g, 60_000);
    const s = await VcfRuntime.open(entry({ index: 'i' }), context({ f: file, i: index }));
    const a = await s.read({ system: 'chr1', start: 0, end: 200_000, level: 1, shapes: ['ranges'] }, signal);
    expect(a.ok && a.tables.ranges?.form === 'list' && a.tables.ranges.start.length).toBeGreaterThan(1000);
    if (a.ok && a.tables.ranges?.form === 'list')
      expect(columnValues(a.tables.ranges.id!).length).toBeGreaterThan(1000);
  });
});
