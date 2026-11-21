/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';

// The genome formats' files: UCSC hg38's cytobands and TP53's region of Ensembl's GFF3 and GENCODE's GTF recorded as
// pages (`fixtures/genomes/sources.json`), and files written by the generators below (a test compares each with its
// generator): BED12 genes on a small generated chromosome, and a misread BED and GFF3 whose bad lines the runtimes say.

/** UCSC hg38's `cytoBandIdeo.txt.gz`, read from its recording (`cytoBandIdeo.txt.gz.pages.json`). */
export const CytobandsFile = 'fixtures/genomes/cytoBandIdeo.txt.gz';

/** Ensembl GRCh38 release 116's GFF3 cut to TP53's region (chromosome `17`), read from its recording. */
export const EnsemblTp53File = 'fixtures/genomes/Homo_sapiens.GRCh38.116.chr.TP53.gff3.gz';
/** GENCODE v50's GTF cut to TP53's region (`chr17`), read from its recording. */
export const GencodeTp53File = 'fixtures/genomes/gencode.v50.annotation.TP53.gtf.gz';
/** TP53's region in GRCh38 (0-based, end exclusive): the gene and its neighbour WRAP53's start. */
export const Tp53Region: readonly [number, number] = [7_660_000, 7_690_000];

/** The examples' files read from recordings (`<uri>.pages.json`): only these are asked for one. */
export const RecordedFiles: readonly string[] = [CytobandsFile, EnsemblTp53File, GencodeTp53File];

/** chr17's length in GRCh38 (the cytobands' last band ends there). */
export const Chr17Length = 83_257_441;

/** BED12 genes on `chrG`, generated. */
export const GeneratedBedFile = 'fixtures/genomes/generated-genes.bed';
export const GeneratedSystem = 'chrG';
export const GeneratedLength = 60_000;

/** A BED with lines the runtime leaves out or reads in part, each said. */
export const MisreadBedFile = 'fixtures/genomes/misread.bed';

/** The generated genes: 24 BED12 lines, exons as blocks, a thick CDS, both strands, sorted. */
export function generatedGenesBed(): string {
  const rng = new PCG(17);
  const lines: string[] = ['track name=genes description="Generated BED12 genes"'];
  let at = 500;
  for (let g = 0; g < 24 && at < GeneratedLength - 3000; g++) {
    const length = 800 + Math.floor(rng.float() * 1600);
    const start = at;
    const end = start + length;
    const exons = 1 + Math.floor(rng.float() * 4);
    const step = Math.floor(length / exons);
    const sizes: number[] = [];
    const starts: number[] = [];
    for (let e = 0; e < exons; e++) {
      const s = e * step;
      starts.push(s);
      sizes.push(e === exons - 1 ? length - s : Math.max(60, Math.floor(step * (0.3 + 0.4 * rng.float()))));
    }
    const strand = rng.float() < 0.5 ? '+' : '-';
    const score = Math.floor(rng.float() * 1000);
    lines.push(
      [
        GeneratedSystem,
        start,
        end,
        `gene${g + 1}`,
        score,
        strand,
        start + Math.floor(sizes[0] / 3),
        end - Math.floor(sizes[exons - 1] / 3),
        '0,0,0',
        exons,
        `${sizes.join(',')},`,
        `${starts.join(',')},`,
      ].join('\t'),
    );
    at = end + 200 + Math.floor(rng.float() * 2000);
  }
  return `${lines.join('\n')}\n`;
}

/** A misread BED6: two good lines and one of each fault the runtime says (a line number each). */
export function misreadBed(): string {
  return [
    '# a BED6 with faults',
    'chrG\t100\t900\tgood1\t500\t+',
    'chrG\t1000\t800\tend-before-start\t0\t+',
    'chrG\tten\t2000\tnot-a-number\t0\t-',
    'chrG\t2100\t2600\tbad-strand\t0\t*',
    'chrG\t3000\t3600\tbad-score\thigh\t-',
    'chrG 4000 4800 spaces 10 +',
    'chrG\t5000\t5900\tgood2\t700\t-\textra',
    '',
  ].join('\n');
}

/** A GFF3 with lines the runtime leaves out or reads flat, each said. */
export const MisreadGffFile = 'fixtures/genomes/misread.gff3';

/** A misread GFF3: a gene with two transcripts read, and one of each fault the runtime says (a line number each). */
export function misreadGff(): string {
  const l = (...c: (string | number)[]) => c.join('\t');
  return [
    '##gff-version 3',
    '##sequence-region chrG 1 7000',
    l('chrG', 'gen', 'gene', 101, 2400, '.', '+', '.', 'ID=gene:A;Name=GOOD'),
    l('chrG', 'gen', 'mRNA', 101, 2400, '.', '+', '.', 'ID=tx:A1;Parent=gene:A;Name=GOOD-201'),
    l('chrG', 'gen', 'exon', 101, 400, '.', '+', '.', 'Parent=tx:A1'),
    l('chrG', 'gen', 'exon', 1801, 2400, '.', '+', '.', 'Parent=tx:A1'),
    l('chrG', 'gen', 'CDS', 201, 400, '.', '+', '0', 'Parent=tx:A1'),
    l('chrG', 'gen', 'CDS', 1801, 2200, '.', '+', '2', 'Parent=tx:A1'),
    l('chrG', 'gen', 'gene', 3000, 2800, '.', '+', '.', 'ID=gene:end-before-start'),
    l('chrG', 'gen', 'gene', 'ten', 3500, '.', '-', '.', 'ID=gene:not-a-number'),
    l('chrG', 'gen', 'gene', 3600, 3900),
    l('chrG', 'gen', 'exon', 4001, 4200, '.', '-', '.', 'Parent=tx:nobody'),
    l('chrG', 'gen', 'mRNA', 4501, 5000, '.', '-', '.', 'ID=tx:loop1;Parent=tx:loop2'),
    l('chrG', 'gen', 'mRNA', 4501, 5000, '.', '-', '.', 'ID=tx:loop2;Parent=tx:loop1'),
    l('chrG', 'gen', 'gene', 5501, 6500, '.', '+', '.', 'ID=gene:B;Name=PROKARYOTE-LIKE'),
    l('chrG', 'gen', 'CDS', 5501, 6500, '.', '+', '0', 'ID=cds:B;Parent=gene:B'),
    '',
  ].join('\n');
}
