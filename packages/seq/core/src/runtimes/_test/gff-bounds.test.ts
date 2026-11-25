/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { columnValues } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { DefaultOpenLimits, type RangesList } from '../../data/protocol.js';
import { assembleTranscripts, GffRuntime, parseGff } from '../gff.js';
import { transcriptTest } from '../gff-assemble.js';
import { GffCauses } from '../gff-issues.js';
import { memoryAssets } from './bed-fixtures.js';

// The GFF runtime's fix round (_spec/README.md "Data", GFF3 and GTF): files that overflowed the stack (a `Parent`
// chain of 20,000 with children before parents, one transcript of 200,000 exons), an encoded comma in `Parent`, a gene
// with an mRNA and a direct exon, a GTF gene's label, `features: genes` one model a gene, a chromosome-long feature
// not making every region read every item, and lines whose separators are counted before they are split.

const signal = new AbortController().signal;
const tsv = (...lines: (string | number)[][]) => table(lines);
/** Lines of cells as text (an array: a long list is never spread). */
const table = (lines: (string | number)[][]) => lines.map((l) => l.join('\t')).join('\n');

function rows(text: string, format: 'gff3' | 'gtf' = 'gff3', causes = new GffCauses()) {
  const f = parseGff(text, format, transcriptTest(), causes, 1e7);
  if (typeof f === 'string') throw new Error(f);
  return f.rows;
}

async function open(text: string, format = 'gff3', features?: string) {
  return GffRuntime.open(
    { address: '/data/s', format, value: text as never, ...(features ? { features } : {}) } as never,
    { assets: memoryAssets({}), signal, limits: DefaultOpenLimits },
  );
}

async function list(text: string, start: number, end: number, features?: string, format = 'gff3') {
  const s = await open(text, format, features);
  const a = await s.read({ system: 'c', start, end, level: 1, shapes: ['ranges'] }, signal);
  if (!a.ok || a.tables.ranges?.form !== 'list') throw new Error('no list');
  return { s, r: a.tables.ranges as RangesList };
}

describe('the GFF runtime’s bounds and fixes', () => {
  it('a `Parent` chain of 20,000, children before parents: no throw, those past 16 flat (`.parent`)', async () => {
    const n = 20_000;
    const lines: (string | number)[][] = [];
    for (let i = n - 1; i >= 0; i--)
      lines.push(['c', '.', 'repeat_region', 1, 100, '.', '+', '.', `ID=r${i}${i ? `;Parent=r${i - 1}` : ''}`]);
    const a = assembleTranscripts(rows(table(lines)));
    expect(a.deep.length).toBe(n - 17);
    const s = await open(table(lines));
    expect(s.issues!.map((i) => i.code)).toContain('seq.source.gff.parent');
  });

  it('one transcript of 200,000 exons: no throw, its gene merges them', () => {
    const lines: (string | number)[][] = [
      ['c', '.', 'gene', 1, 2_000_000, '.', '+', '.', 'ID=g'],
      ['c', '.', 'mRNA', 1, 2_000_000, '.', '+', '.', 'ID=t;Parent=g'],
    ];
    for (let e = 0; e < 200_000; e++) lines.push(['c', '.', 'exon', e * 10 + 1, e * 10 + 5, '.', '+', '.', 'Parent=t']);
    const a = assembleTranscripts(rows(table(lines)));
    expect(a.genes[0].pieces.length).toBe(400_000);
  });

  it('an encoded comma in `Parent` and `Alias` is part of the value, not a separator', async () => {
    const text = tsv(
      ['c', '.', 'mRNA', 1, 100, '.', '+', '.', 'ID=t%2C1;Alias=A%2CB,C'],
      ['c', '.', 'exon', 1, 50, '.', '+', '.', 'Parent=t%2C1'],
    );
    const r = rows(text);
    expect(r[1].parents).toEqual(['t,1']);
    const a = assembleTranscripts(r);
    expect(a.transcripts.map((t) => [t.id, t.pieces])).toEqual([['t,1', [0, 50]]]);
    const s = await open(text);
    expect((await s.findNames!('a,b', { max: 5 }, signal)).map((h) => h.field)).toEqual(['Alias']);
    expect((await s.findNames!('c', { max: 5 }, signal)).map((h) => h.field)).toEqual(['Alias']);
  });

  it('a gene with an mRNA and a direct exon: the mRNA keeps its gene', () => {
    const a = assembleTranscripts(
      rows(
        tsv(
          ['c', '.', 'gene', 1, 900, '.', '+', '.', 'ID=g;Name=G'],
          ['c', '.', 'mRNA', 1, 900, '.', '+', '.', 'ID=t;Parent=g'],
          ['c', '.', 'exon', 1, 100, '.', '+', '.', 'Parent=t'],
          ['c', '.', 'exon', 801, 900, '.', '+', '.', 'Parent=g'],
        ),
      ),
    );
    expect(a.transcripts.map((t) => [t.id, t.gene, t.geneName])).toEqual([['t', 'g', 'G']]);
    expect(a.genes.map((g) => g.id)).toEqual(['g']);
  });

  it('a GTF gene found by `gene_id` is labelled by its `gene_name`, not a transcript’s', async () => {
    const text = tsv(
      ['c', '.', 'gene', 1, 900, '.', '+', '.', 'gene_id "G1"; gene_name "ALPHA";'],
      [
        'c',
        '.',
        'transcript',
        1,
        900,
        '.',
        '+',
        '.',
        'gene_id "G1"; transcript_id "T1"; gene_name "ALPHA"; transcript_name "ALPHA-201";',
      ],
      ['c', '.', 'exon', 1, 100, '.', '+', '.', 'gene_id "G1"; transcript_id "T1";'],
    );
    const s = await open(text, 'gtf');
    const hits = await s.findNames!('g1', { max: 5 }, signal);
    expect(hits.map((h) => [h.field, h.label])).toEqual([['gene_id', 'ALPHA']]);
  });

  it('`features: genes` is one model a gene: no other feature an item', async () => {
    const text = tsv(
      ['c', '.', 'gene', 1, 900, '.', '+', '.', 'ID=g'],
      ['c', '.', 'mRNA', 1, 900, '.', '+', '.', 'ID=t;Parent=g'],
      ['c', '.', 'exon', 1, 100, '.', '+', '.', 'Parent=t'],
      ['c', '.', 'biological_region', 200, 300, '.', '+', '.', 'logic_name=x'],
    );
    const { r } = await list(text, 0, 1000, 'genes');
    expect(new Set(columnValues(r.id!))).toEqual(new Set(['g']));
  });

  it('`features: rows` collects an `ID`’s rows once (linear): 20,000 CDS rows of one ID', async () => {
    const lines: (string | number)[][] = [];
    for (let e = 0; e < 20_000; e++) lines.push(['c', '.', 'CDS', e * 10 + 1, e * 10 + 5, '.', '+', '0', 'ID=cds1']);
    const { r } = await list(table(lines), 0, 50, 'rows');
    expect(columnValues(r.id!)).toEqual(['cds1', 'cds1', 'cds1', 'cds1', 'cds1']);
    expect(r.itemEnd![0]).toBe(199_995);
  });

  it('a chromosome-long feature: a region far from the start reads only what meets it', async () => {
    const lines: (string | number)[][] = [['c', '.', 'repeat_region', 1, 10_000_000, '.', '+', '.', 'ID=long']];
    for (let i = 0; i < 2000; i++)
      lines.push(['c', '.', 'repeat_region', i * 5000 + 1, i * 5000 + 100, '.', '+', '.', `ID=r${i}`]);
    const { r } = await list(table(lines), 5_000_000, 5_000_200);
    expect(new Set(columnValues(r.id!))).toEqual(new Set(['long', 'r1000']));
  });

  it('a line of a million tabs or semicolons is refused or bounded without splitting it whole', () => {
    const causes = new GffCauses();
    const tabs = `c\t.\tgene\t1\t10\t.\t+\t.\tID=a${'\t'.repeat(1_000_000)}`;
    const semis = `c\t.\tgene\t1\t10\t.\t+\t.\tID=b${Array.from({ length: 300 }, (_, i) => `;k${i}=v`).join('')}`;
    const r = rows(`${tabs}\n${semis}`, 'gff3', causes);
    expect(r.map((x) => x.id)).toEqual(['b']);
    expect(Object.keys(r[0].attributes).length).toBeLessThanOrEqual(256);
    expect(causes.issues('/s').map((i) => i.code)).toEqual(['seq.source.gff.line', 'seq.source.gff.attributes']);
  });
});
