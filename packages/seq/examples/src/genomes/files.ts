/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';

// The genome formats' files: UCSC hg38's cytobands, TP53's region of Ensembl's GFF3 and GENCODE's GTF, and ClinVar's
// sites file (TP53's pages) with its index recorded as pages (`fixtures/genomes/sources.json`), and files written by the
// generators below (a test compares each with its generator): BED12 genes and variants on small generated chromosomes,
// and a misread BED, GFF3 and VCF whose bad lines the runtimes say.

/** UCSC hg38's `cytoBandIdeo.txt.gz`, read from its recording (`cytoBandIdeo.txt.gz.pages.json`). */
export const CytobandsFile = 'fixtures/genomes/cytoBandIdeo.txt.gz';

/** Ensembl GRCh38 release 116's GFF3 cut to TP53's region (chromosome `17`), read from its recording. */
export const EnsemblTp53File = 'fixtures/genomes/Homo_sapiens.GRCh38.116.chr.TP53.gff3.gz';
/** GENCODE v50's GTF cut to TP53's region (`chr17`), read from its recording. */
export const GencodeTp53File = 'fixtures/genomes/gencode.v50.annotation.TP53.gtf.gz';
/** TP53's region in GRCh38 (0-based, end exclusive): the gene and its neighbour WRAP53's start. */
export const Tp53Region: readonly [number, number] = [7_660_000, 7_690_000];

/** NCBI ClinVar's GRCh38 sites of 2024-01-07 (BGZF), read by region from its recording: TP53's pages only. */
export const ClinvarFile = 'fixtures/genomes/clinvar_20240107.vcf.gz';
/** ClinVar's tabix index, named by the source (never `<uri>.tbi`), recorded whole. */
export const ClinvarIndexFile = 'fixtures/genomes/clinvar_20240107.vcf.gz.tbi';

/** The examples' files read from recordings (`<uri>.pages.json`): only these are asked for one. */
export const RecordedFiles: readonly string[] = [
  CytobandsFile,
  EnsemblTp53File,
  GencodeTp53File,
  ClinvarFile,
  ClinvarIndexFile,
];

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

/** Variants on `chrG` (20 kb), generated: SNVs, a few multi-allelic, indels, and two long deletions over SNVs. */
export const GeneratedVcfFile = 'fixtures/genomes/generated-variants.vcf';
export const GeneratedVcfLength = 20_000;

/** The generated variants: a VCF 4.2 with its contig's length, sorted, every category the runtime reads. */
export function generatedVariantsVcf(): string {
  const rng = new PCG(23);
  const base = () => 'ACGT'[Math.floor(rng.float() * 4)];
  const other = (b: string) => 'ACGT'.replace(b, '')[Math.floor(rng.float() * 3)];
  const rows: [number, string][] = [];
  const line = (pos: number, id: string, ref: string, alt: string, info: string) =>
    rows.push([pos, [GeneratedSystem, pos, id, ref, alt, '.', 'PASS', info].join('\t')]);
  for (let pos = 150, k = 1; pos < GeneratedVcfLength - 200; pos += 120 + Math.floor(rng.float() * 500), k++) {
    const r = rng.float();
    const ref = base();
    const af = (x: number) => (Math.floor(x * 1000) / 1000).toFixed(3);
    if (r < 0.6) line(pos, `snv${k}`, ref, other(ref), `AF=${af(rng.float() / 2)}`);
    else if (r < 0.7) {
      const a = other(ref);
      const b = 'ACGT'.replace(ref, '').replace(a, '')[0];
      line(pos, `multi${k}`, ref, `${a},${b}`, `AF=${af(rng.float() / 4)},${af(rng.float() / 4)}`);
    } else if (r < 0.8) line(pos, `ins${k}`, ref, ref + base() + base() + base(), 'AF=0.010');
    else if (r < 0.9) line(pos, `del${k}`, ref + base() + base() + base() + base(), ref, 'AF=0.020');
    else line(pos, `mnv${k}`, ref + base(), other(ref) + base(), 'AF=0.005');
  }
  // two long deletions, each one item across the SNVs it covers: a symbolic one to its END, a sequence-resolved one
  line(4001, 'sv1', 'N', '<DEL>', 'SVTYPE=DEL;END=9000;AF=0.100');
  let long = '';
  for (let i = 0; i < 600; i++) long += base();
  line(12_001, 'sv2', long, long[0], 'AF=0.050');
  rows.sort((a, b) => a[0] - b[0]);
  return `${[
    '##fileformat=VCFv4.2',
    `##contig=<ID=${GeneratedSystem},length=${GeneratedVcfLength}>`,
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the structural variant">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End of the variant">',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    ...rows.map((r) => r[1]),
  ].join('\n')}\n`;
}

/** A VCF with lines the runtime leaves out or reads in part, each said. */
export const MisreadVcfFile = 'fixtures/genomes/misread.vcf';

/** A misread VCF: good lines and one of each fault the runtime says (a line number each). */
export function misreadVcf(): string {
  const l = (...c: (string | number)[]) => c.join('\t');
  return [
    '##fileformat=VCFv4.2',
    '##contig=<ID=chrG,length=7000>',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##INFO=<ID=BROKEN>',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    l('chrG', 500, 'good1', 'A', 'G', 50, 'PASS', 'DP=20'),
    l('chrG', 'ten', 'not-a-number', 'A', 'G', 50, 'PASS', 'DP=20'),
    l('chrG', 1500, 'no-ref', '', 'G', 50, 'PASS', 'DP=20'),
    l('chrG', 2000, 'short-line', 'A'),
    l('chrG', 2500, 'bad-depth', 'C', 'T', 'high', 'PASS', 'DP=many;XX=1'),
    l('chrG', 3000, 'star', 'G', 'C,*', 30, 'PASS', 'DP=8'),
    l('chrG', 6990, 'past-end', 'ACGTACGTACGTACGTACGT', 'A', 40, 'PASS', 'DP=5'),
    l('chrG', 7500, 'outside', 'A', 'T', 40, 'PASS', 'DP=5'),
    l('chrG', 5000, 'good2', 'T', 'TA', 60, 'PASS', 'DP=30'),
    '',
  ].join('\n');
}
