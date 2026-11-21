/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { gzipSync } from 'node:zlib';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { columnValues } from '@molstar/core/repr/tree';
import { bgzip } from '@molstar/test-harness/genome/bgzf-writer';
import { describe, expect, it } from 'vitest';
import { sourceSuite } from '../../data/_test/source-suite.js';
import { delayedRuntime } from '../../data/async-source.js';
import { findNames } from '../../data/names.js';
import {
  DefaultOpenLimits,
  type OpenContext,
  type OpenedSource,
  type RangesList,
  type RegionRequest,
  type SourceEntry,
} from '../../data/protocol.js';
import { SourceOpenError } from '../../data/seq-runtime.js';
import { checkAnswer } from '../../data/tables.js';
import { assembleTranscripts, GffRuntime, gffRuntime, parseGff } from '../gff.js';
import { GffCauses } from '../gff-issues.js';
import { transcriptTest } from '../gff-assemble.js';
import { memoryAssets } from './bed-fixtures.js';
import { EnsemblTp53, GencodeTp53, generateGff3, generateGtf, recordedSlice } from './gff-fixtures.js';

// GFF3 and GTF (_spec/README.md "Data", GFF3 and GTF): the source suite on generated files and on the recorded real
// slices (TP53 from Ensembl's GFF3 and GENCODE's GTF); `assembleTranscripts` alone (three levels, a shared exon,
// CDS-only, GTF with and without transcript lines, a cycle, a chain past 16, a prokaryote's gene, a type outside the
// list); the `features` modes; span rows; lengths; the name index; fuzzed and truncated files without a throw.

const signal = new AbortController().signal;
const names = ['chr1', 'chr2'];

function context(files: Record<string, Uint8Array | string>, limits = DefaultOpenLimits): OpenContext {
  return { assets: memoryAssets(files), signal, limits };
}

const entry = (format: string, more: Partial<SourceEntry> & { features?: string } = {}): SourceEntry => ({
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

/** Rows of a GFF3 text (`##gff-version 3` implied), parsed as the runtime does. */
function rows(text: string, format: 'gff3' | 'gtf' = 'gff3') {
  const f = parseGff(text, format, transcriptTest(), new GffCauses(), 1e6);
  if (typeof f === 'string') throw new Error(f);
  return f.rows;
}

const gff = (...lines: string[][]) => lines.map((l) => l.join('\t')).join('\n');

async function read(s: OpenedSource, system: string, start: number, end: number): Promise<RangesList> {
  const request: RegionRequest = { system, start, end, level: 1, shapes: ['ranges'] };
  const a = await s.read(request, signal);
  expect(checkAnswer(a, request)).toBeUndefined();
  if (!a.ok || a.tables.ranges?.form !== 'list') throw new Error('no list');
  return a.tables.ranges;
}

const field = (r: RangesList, name: string) => columnValues(r.fields!.find((f) => f.name === name)!.column);

describe('the source suite', () => {
  const g3 = generateGff3(60, names, 3);
  const gtf = generateGtf(60, names, 3);

  it('GFF3 read whole: plain, gzip and BGZF by their magic; every features mode; the synthetic twin', async () => {
    const text = new TextEncoder().encode(g3.text);
    for (const bytes of [text, new Uint8Array(gzipSync(text)), bgzip(text, 3000).bytes])
      await sourceSuite(GffRuntime, entry('gff3'), context({ f: bytes }));
    for (const features of ['genes', 'rows'])
      await sourceSuite(GffRuntime, entry('gff3', { features }), context({ f: g3.text }));
    await sourceSuite(delayedRuntime(GffRuntime, { delayMs: 1 }), entry('gff3'), context({ f: g3.text }));
  });

  it('GTF read whole, with and without transcript lines', async () => {
    await sourceSuite(GffRuntime, entry('gtf'), context({ f: gtf.text }));
    await sourceSuite(GffRuntime, entry('gtf'), context({ f: generateGtf(60, names, 3, false).text }));
  });

  it('the recorded real slices: TP53 from Ensembl’s GFF3 and GENCODE’s GTF, with no issue but the description', async () => {
    for (const [format, name, system] of [
      ['gff3', EnsemblTp53, '17'],
      ['gtf', GencodeTp53, 'chr17'],
    ] as const) {
      const bytes = recordedSlice(name);
      await sourceSuite(GffRuntime, entry(format), context({ f: bytes }));
      const s = await GffRuntime.open(entry(format), context({ f: bytes }));
      expect(s.issues!.map((i) => i.code).filter((c) => c !== 'seq.source.gff.region')).toEqual([]);
      const r = await read(s, system, 7_660_000, 7_700_000);
      const labels = new Set(field(r, 'label'));
      // TP53's canonical transcript, by its name; the transcripts' gene named
      expect(labels.has('TP53-201')).toBe(true);
      expect(new Set(field(r, 'gene_name')).has('TP53')).toBe(true);
      // a coding transcript's thick part inside its extent
      const i = field(r, 'label').indexOf('TP53-201');
      const ts = field(r, 'thick_start')[i] as number;
      expect(ts).toBeGreaterThanOrEqual(r.itemStart![i]);
      expect(field(r, 'thick_end')[i] as number).toBeLessThanOrEqual(r.itemEnd![i]);
      expect(r.strand![i]).toBe(-1);
      s.dispose();
    }
  });

  it('the recorded GFF3 slice: chromosome 17’s length from `##sequence-region` (EC-188), its chromosome line no item', async () => {
    const s = await GffRuntime.open(entry('gff3'), context({ f: recordedSlice(EnsemblTp53) }));
    expect(s.coverage.systems).toEqual([{ name: '17', length: 83_257_441 }]);
    expect(s.issues!.map((i) => i.code)).toEqual(['seq.source.gff.region']);
    const r = await read(s, '17', 0, 83_257_441 > 1 << 30 ? 0 : 10_000_000);
    expect(field(r, 'type')).not.toContain('chromosome');
  });
});

describe('assembleTranscripts (EC-227)', () => {
  it('three levels: gene → mRNA → exon, CDS and UTRs; the thick part the CDS hull, the phase of the first CDS', () => {
    const a = assembleTranscripts(
      rows(
        gff(
          ['c', '.', 'gene', '100', '900', '.', '-', '.', 'ID=g1;Name=G'],
          ['c', '.', 'mRNA', '100', '900', '.', '-', '.', 'ID=t1;Parent=g1;Name=G-1'],
          ['c', '.', 'exon', '100', '300', '.', '-', '.', 'Parent=t1'],
          ['c', '.', 'exon', '600', '900', '.', '-', '.', 'Parent=t1'],
          ['c', '.', 'CDS', '200', '300', '.', '-', '2', 'Parent=t1'],
          ['c', '.', 'CDS', '600', '800', '.', '-', '0', 'Parent=t1'],
          ['c', '.', 'three_prime_UTR', '100', '199', '.', '-', '.', 'Parent=t1'],
          ['c', '.', 'five_prime_UTR', '801', '900', '.', '-', '.', 'Parent=t1'],
        ),
      ),
    );
    expect(a.transcripts).toHaveLength(1);
    const t = a.transcripts[0];
    expect(t).toMatchObject({ id: 't1', start: 99, end: 900, strand: -1, gene: 'g1', geneName: 'G' });
    expect(t.pieces).toEqual([99, 300, 599, 900]);
    expect([t.thickStart, t.thickEnd]).toEqual([199, 800]);
    // reversed: the first CDS is the one ending last
    expect(t.cdsPhase).toBe(0);
    expect(a.genes.map((g) => g.id)).toEqual(['g1']);
    expect(a.features).toEqual([]);
  });

  it('a shared exon (`Parent=t1,t2`) is a piece of each', () => {
    const a = assembleTranscripts(
      rows(
        gff(
          ['c', '.', 'mRNA', '1', '500', '.', '+', '.', 'ID=t1'],
          ['c', '.', 'mRNA', '1', '800', '.', '+', '.', 'ID=t2'],
          ['c', '.', 'exon', '1', '100', '.', '+', '.', 'Parent=t1,t2'],
          ['c', '.', 'exon', '400', '500', '.', '+', '.', 'Parent=t1'],
          ['c', '.', 'exon', '700', '800', '.', '+', '.', 'Parent=t2'],
        ),
      ),
    );
    expect(a.transcripts.map((t) => [t.id, t.pieces])).toEqual([
      ['t1', [0, 100, 399, 500]],
      ['t2', [0, 100, 699, 800]],
    ]);
  });

  it('CDS only (no exon lines): the pieces its CDS and UTRs merged', () => {
    const a = assembleTranscripts(
      rows(
        gff(
          ['c', '.', 'mRNA', '1', '500', '.', '+', '.', 'ID=t1'],
          ['c', '.', 'five_prime_UTR', '1', '49', '.', '+', '.', 'Parent=t1'],
          ['c', '.', 'CDS', '50', '200', '.', '+', '0', 'Parent=t1'],
          ['c', '.', 'CDS', '300', '500', '.', '+', '1', 'Parent=t1'],
        ),
      ),
    );
    expect(a.transcripts[0].pieces).toEqual([0, 200, 299, 500]);
    expect(a.transcripts[0].cdsPhase).toBe(0);
  });

  it('GTF with transcript lines: the extent the line’s; without: the hull of its rows, made from them', () => {
    const lines = (withTx: boolean) =>
      gff(
        ['c', '.', 'gene', '1', '1000', '.', '+', '.', 'gene_id "G"; gene_name "GN";'],
        ...(withTx ? [['c', '.', 'transcript', '1', '1000', '.', '+', '.', 'gene_id "G"; transcript_id "T";']] : []),
        ['c', '.', 'exon', '10', '100', '.', '+', '.', 'gene_id "G"; transcript_id "T"; tag "a"; tag "b";'],
        ['c', '.', 'exon', '900', '950', '.', '+', '.', 'gene_id "G"; transcript_id "T";'],
        ['c', '.', 'CDS', '50', '100', '.', '+', '0', 'gene_id "G"; transcript_id "T";'],
        ['c', '.', 'start_codon', '50', '52', '.', '+', '0', 'gene_id "G"; transcript_id "T";'],
        ['c', '.', 'stop_codon', '900', '902', '.', '+', '0', 'gene_id "G"; transcript_id "T";'],
      );
    const w = assembleTranscripts(rows(lines(true), 'gtf'));
    expect(w.transcripts[0]).toMatchObject({ id: 'T', start: 0, end: 1000, gene: 'G', geneName: 'GN' });
    expect(w.synthetic).toBe(0);
    const wo = assembleTranscripts(rows(lines(false), 'gtf'));
    expect(wo.transcripts[0]).toMatchObject({ id: 'T', start: 9, end: 950, gene: 'G', synthetic: true });
    expect(wo.synthetic).toBe(1);
    // GTF's stop codon outside the CDS widens the thick part; a repeated key joined
    expect([wo.transcripts[0].thickStart, wo.transcripts[0].thickEnd]).toEqual([49, 902]);
    expect(rows(lines(false), 'gtf')[1].attributes.tag).toBe('a,b');
  });

  it('a cycle and a chain deeper than 16 leave their rows items of their own', () => {
    const cyc = assembleTranscripts(
      rows(
        gff(
          ['c', '.', 'mRNA', '1', '100', '.', '+', '.', 'ID=a;Parent=b'],
          ['c', '.', 'region', '1', '100', '.', '+', '.', 'ID=b;Parent=a'],
          ['c', '.', 'exon', '1', '100', '.', '+', '.', 'Parent=a'],
        ),
      ),
    );
    expect(cyc.cycles).toEqual([1, 2]);
    expect(cyc.transcripts).toEqual([]);
    expect(cyc.features.map((f) => f.id)).toEqual(['a', 'b', 'exon:1-100']);
    const chain = Array.from({ length: 20 }, (_, i) => [
      'c',
      '.',
      'region',
      '1',
      '100',
      '.',
      '+',
      '.',
      `ID=r${i}${i ? `;Parent=r${i - 1}` : ''}`,
    ]);
    const deep = assembleTranscripts(rows(gff(...chain)));
    expect(deep.deep.length).toBe(3);
    expect(deep.features).toHaveLength(20);
  });

  it('a prokaryote’s gene → CDS makes the gene a transcript-level item (counted by type); a type outside the list too', () => {
    const a = assembleTranscripts(
      rows(
        gff(
          ['c', '.', 'gene', '1', '300', '.', '+', '.', 'ID=g;Name=dnaA'],
          ['c', '.', 'CDS', '1', '300', '.', '+', '0', 'ID=cds;Parent=g'],
          ['c', '.', 'mystery_RNA', '400', '600', '.', '-', '.', 'ID=m'],
          ['c', '.', 'exon', '400', '450', '.', '-', '.', 'Parent=m'],
          ['c', '.', 'exon', '550', '600', '.', '-', '.', 'Parent=m'],
          ['c', '.', 'tmRNA_transcript', '700', '800', '.', '+', '.', 'ID=x'],
        ),
      ),
    );
    expect(a.transcripts.map((t) => [t.id, t.pieces])).toEqual([
      ['g', [0, 300]],
      ['m', [399, 450, 549, 600]],
      ['x', [699, 800]],
    ]);
    expect(a.transcripts[0]).toMatchObject({ thickStart: 0, thickEnd: 300 });
    expect([...a.promoted]).toEqual([
      ['gene', 1],
      ['mystery_RNA', 1],
    ]);
  });

  it('attributes past 256 dropped, said', async () => {
    const attrs = Array.from({ length: 300 }, (_, i) => `k${i}=v`).join(';');
    const s = await GffRuntime.open(
      entry('gff3', {
        uri: undefined,
        value: gff(['c', '.', 'gene', '1', '10', '.', '+', '.', `ID=g;${attrs}`]) as never,
      }),
      context({}),
    );
    expect(s.issues!.map((i) => i.code)).toContain('seq.source.gff.attributes');
    expect(s.coverage.shapes[0].fields!.length).toBeLessThan(270);
  });
});

describe('the runtime', () => {
  const text = gff(
    ['c', '.', 'gene', '100', '900', '.', '+', '.', 'ID=g1;Name=ALPHA;Alias=A1,AL'],
    ['c', '.', 'mRNA', '100', '900', '.', '+', '.', 'ID=t1;Parent=g1;Name=ALPHA-201'],
    ['c', '.', 'mRNA', '100', '700', '.', '+', '.', 'ID=t2;Parent=g1;Name=ALPHA-202'],
    ['c', '.', 'exon', '100', '200', '.', '+', '.', 'Parent=t1,t2'],
    ['c', '.', 'exon', '800', '900', '.', '+', '.', 'Parent=t1'],
    ['c', '.', 'exon', '600', '700', '.', '+', '.', 'Parent=t2'],
    ['c', '.', 'CDS', '150', '200', '.', '+', '0', 'Parent=t1'],
    ['c', '.', 'CDS', '800', '850', '.', '+', '1', 'Parent=t1'],
    ['c', '.', 'repeat_region', '950', '990', '.', '.', '.', 'Name=rep'],
  );
  const open = (more: Partial<SourceEntry> & { features?: string } = {}) =>
    GffRuntime.open(entry('gff3', { uri: undefined, value: text as never, ...more }), context({}));

  it('transcripts (the default): a transcript an item, its pieces with its extent; other features items of their own', async () => {
    const s = await open();
    const r = await read(s, 'c', 0, 1000);
    expect(columnValues(r.id!)).toEqual(['t1', 't2', 't2', 't1', 'repeat_region:950-990']);
    expect(Array.from(r.itemStart!)).toEqual([99, 99, 99, 99, 949]);
    expect(Array.from(r.itemEnd!)).toEqual([900, 700, 700, 900, 990]);
    expect(field(r, 'gene')).toEqual(['g1', 'g1', 'g1', 'g1', null]);
    expect(field(r, 'label')).toEqual(['ALPHA-201', 'ALPHA-202', 'ALPHA-202', 'ALPHA-201', 'rep']);
    expect(field(r, 'cds_phase')).toEqual([0, -1, -1, 0, -1]);
    expect(s.coverage.names).toEqual(['Name', 'ID', 'Alias']);
  });

  it('a region inside an intron gives each transcript there a span row (EC-229)', async () => {
    const s = await open();
    const r = await read(s, 'c', 300, 400);
    expect(columnValues(r.id!)).toEqual(['t1', 't2']);
    expect(Array.from(r.start)).toEqual([300, 300]);
    expect(Array.from(r.end)).toEqual([300, 300]);
    expect(Array.from(r.itemEnd!)).toEqual([900, 700]);
  });

  it('genes: one item a gene, its transcripts’ exons merged, its CDS hull; rows: the file’s lines', async () => {
    const g = await read(await open({ features: 'genes' } as never), 'c', 0, 1000);
    expect(columnValues(g.id!)).toEqual(['g1', 'g1', 'g1', 'repeat_region:950-990']);
    expect(Array.from(g.start)).toEqual([99, 599, 799, 949]);
    expect(field(g, 'thick_start')[0]).toBe(149);
    const r = await read(await open({ features: 'rows' } as never), 'c', 0, 1000);
    expect(r.start.length).toBe(9);
    expect(field(r, 'type')).toContain('CDS');
    expect(field(r, 'phase')).toContain(1);
  });

  it('the fixed conventions said once if written; an `index` said, the file read whole', async () => {
    const s = await open({ base: 0, ends: 'exclusive', index: 'x.tbi' });
    expect(s.issues!.map((i) => i.code).sort()).toEqual(['seq.source.gff.index', 'seq.source.gff.transform']);
  });

  it('`##FASTA` and what follows ignored, counted; `##sequence-region` a length (EC-188); a whole-sequence region no item', async () => {
    const t = [
      '##gff-version 3',
      '##sequence-region ctg 1 5000',
      gff(
        ['ctg', '.', 'region', '1', '5000', '.', '+', '.', 'ID=ctg'],
        ['ctg', '.', 'gene', '10', '90', '.', '+', '.', 'ID=g'],
      ),
      '##FASTA',
      '>ctg',
      'ACGT'.repeat(100),
    ].join('\n');
    const s = await GffRuntime.open(entry('gff3', { uri: undefined, value: t as never }), context({}));
    expect(s.coverage.systems).toEqual([{ name: 'ctg', length: 5000 }]);
    expect(s.issues!.map((i) => i.code).sort()).toEqual(['seq.source.gff.fasta', 'seq.source.gff.region']);
    const r = await read(s, 'ctg', 0, 5000);
    expect(columnValues(r.id!)).toEqual(['g']);
  });

  it('lines that do not parse are left out with their numbers; a missing parent and a cycle said', async () => {
    const t = gff(
      ['c', '.', 'gene', '10', '5', '.', '+', '.', 'ID=bad'],
      ['c', '.', 'gene', 'x', '50', '.', '+', '.', 'ID=bad2'],
      ['c', '.', 'gene'],
      ['c', '.', 'exon', '10', '50', '.', '+', '.', 'Parent=nobody'],
      ['c', '.', 'mRNA', '1', '100', '.', '+', '.', 'ID=a;Parent=b'],
      ['c', '.', 'mRNA', '1', '100', '.', '+', '.', 'ID=b;Parent=a'],
    );
    const s = await GffRuntime.open(entry('gff3', { uri: undefined, value: t as never }), context({}));
    const issues = s.issues!.map((i) => `${i.code} ${i.message}`);
    expect(issues.some((m) => m.startsWith('seq.source.gff.line An end before its start: 1 (line 1)'))).toBe(true);
    expect(issues.filter((m) => m.startsWith('seq.source.gff.line')).length).toBe(3);
    expect(issues.some((m) => m.includes('A `Parent` no line has: 1 (line 4)'))).toBe(true);
    expect(issues.some((m) => m.includes('A cycle'))).toBe(true);
  });

  it('a GTF with an `index` and no transcript lines says its extents are partial by region', async () => {
    const t = generateGtf(4, ['c'], 1, false).text;
    const s = await GffRuntime.open(entry('gtf', { uri: undefined, value: t as never, index: 'i' }), context({}));
    expect(s.issues!.map((i) => i.code)).toContain('seq.source.gtf.extent');
  });

  it('`transcriptTypes` a runtime option', async () => {
    const t = gff(
      ['c', '.', 'odd', '1', '100', '.', '+', '.', 'ID=o'],
      ['c', '.', 'exon', '1', '50', '.', '+', '.', 'Parent=o'],
    );
    const s = await gffRuntime({ transcriptTypes: ['odd'] }).open(
      entry('gff3', { uri: undefined, value: t as never }),
      context({}),
    );
    expect(s.issues).toEqual([]);
  });

  it('bounds: rows past `limits.maxRows` and bytes past `limits.maxBytes` refuse the source, hinting', async () => {
    const g = generateGff3(40, names, 2);
    expect(
      await opened(GffRuntime.open(entry('gff3'), context({ f: g.text }, { maxRows: 100, maxBytes: 1e8 }))),
    ).toEqual(['seq.source.too_large']);
    expect(
      await opened(GffRuntime.open(entry('gff3'), context({ f: g.text }, { maxRows: 1e6, maxBytes: 1000 }))),
    ).toEqual(['seq.source.too_large']);
    expect(await opened(GffRuntime.open(entry('gff3', { uri: 'missing' }), context({})))).toEqual([
      'seq.source.unreadable',
    ]);
  });
});

describe('the name index (EC-228)', () => {
  it('GFF3: `Name`, `ID` and `Alias`, case-insensitive, exact first then by prefix, at most `max`; an id only for an item served', async () => {
    const s = await GffRuntime.open(entry('gff3'), context({ f: recordedSlice(EnsemblTp53) }));
    const exact = await s.findNames!('tp53', { max: 20 }, signal);
    // the gene is no item of `transcripts`: found by its range
    expect(exact[0]).toMatchObject({ system: '17', label: 'TP53', field: 'Name' });
    expect(exact[0].id).toBeUndefined();
    const prefix = await s.findNames!('TP53-2', { max: 5, prefix: true }, signal);
    expect(prefix).toHaveLength(5);
    expect(prefix.every((h) => h.label.startsWith('TP53-2') && h.id?.startsWith('transcript:'))).toBe(true);
    const byId = await s.findNames!('gene:ENSG00000141510', { max: 20 }, signal);
    expect(byId[0]).toMatchObject({ field: 'ID', label: 'TP53' });
    const alias = await GffRuntime.open(
      entry('gff3', {
        uri: undefined,
        value: gff(['c', '.', 'gene', '5', '50', '.', '+', '.', 'ID=g;Name=X;Alias=Old1,Old2']) as never,
      }),
      context({}),
    );
    expect(await alias.findNames!('old2', { max: 20 }, signal)).toEqual([
      { system: 'c', start: 4, end: 50, id: 'g', label: 'X', field: 'Alias' },
    ]);
    // `features: genes` serves the gene: its hit has the id
    const genes = await GffRuntime.open(
      entry('gff3', { features: 'genes' } as never),
      context({ f: recordedSlice(EnsemblTp53) }),
    );
    expect((await genes.findNames!('TP53', { max: 1 }, signal))[0].id).toBe('gene:ENSG00000141510');
  });

  it('GTF: `gene_name`, `gene_id`, `transcript_name`, `transcript_id`; merged across sources by `findNames`', async () => {
    const gtf = await GffRuntime.open(entry('gtf'), context({ f: recordedSlice(GencodeTp53) }));
    expect(gtf.coverage.names).toEqual(['gene_name', 'gene_id', 'transcript_name', 'transcript_id']);
    const t = await gtf.findNames!('tp53-201', { max: 20 }, signal);
    expect(t).toHaveLength(1);
    expect(t[0]).toMatchObject({ field: 'transcript_name', label: 'TP53-201', system: 'chr17' });
    expect(t[0].id).toMatch(/^ENST/);
    const g3 = await GffRuntime.open(entry('gff3'), context({ f: recordedSlice(EnsemblTp53) }));
    const merged = await findNames(
      [
        { address: '/data/gff', source: g3 },
        { address: '/data/gtf', source: gtf },
      ],
      'WRAP53',
      { max: 20 },
      signal,
    );
    expect(merged.map((h) => h.source)).toContain('/data/gff');
    expect(merged.map((h) => h.source)).toContain('/data/gtf');
    expect(await gtf.findNames!('nothing-like-this', { max: 20, prefix: true }, signal)).toEqual([]);
  });
});

describe('fuzzed and truncated files (EC-40, EC-92)', () => {
  it('read without a throw: every outcome a source or a refusal', async () => {
    const rng = new PCG(5);
    const base = new TextEncoder().encode(generateGff3(20, names, 9).text);
    const gtf = new TextEncoder().encode(generateGtf(20, names, 9).text);
    for (let k = 0; k < 60; k++) {
      const src = k % 2 ? gtf : base;
      const b = src.slice(0, Math.floor(rng.float() * src.length));
      for (let m = 0; m < 20; m++) b[Math.floor(rng.float() * b.length)] = Math.floor(rng.float() * 256);
      const got = await opened(GffRuntime.open(entry(k % 2 ? 'gtf' : 'gff3'), context({ f: b })));
      expect(got.length).toBeGreaterThan(0);
    }
    // a gzip cut short is refused (its size read from the cut: too large, or malformed), never a throw
    const gz = new Uint8Array(gzipSync(base));
    const cut = await opened(GffRuntime.open(entry('gff3'), context({ f: gz.subarray(0, gz.length - 9) })));
    expect(['seq.source.malformed', 'seq.source.too_large']).toContain(cut[0]);
  });
});

describe('size', () => {
  // the quick representative of 100,000 lines under 500 ms (`gff-time.test.ts`, heavy): no clock here
  it('10,000 GFF3 lines assembled, every transcript once', () => {
    const g = generateGff3(700, names, 4, 30_000_000);
    expect(g.lines).toBeGreaterThan(10_000);
    const a = assembleTranscripts(rows(g.text));
    expect(a.transcripts.length).toBe(g.transcripts);
    expect(a.genes.length).toBe(g.genes);
  });
});
