/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ExampleView } from '../views.js';
import * as genomeBed from './bed.js';
import * as genomeGff from './gff.js';
import * as genomeVcf from './vcf.js';

/** Views made for the genome formats' stories (`Features/Genome formats`): files read by the opt-in genome runtimes. */
export const GenomeFormatsViews: readonly ExampleView[] = [
  {
    dataset: 'genomes',
    id: 'bed-cytobands',
    title: 'hg38 cytobands from BED',
    summary:
      'chr17’s stained bands from UCSC’s cytoBandIdeo (hg38, the recorded file, gzip): a `bed 4+1` file whose fifth column is named by its number (`stain: 4`), drawn as bands and as named ranges, and as the ideogram above showing the window. Zoom in to a band; hover one for its name and stain.',
    build: genomeBed.cytobands,
    registry: genomeBed.registry,
  },
  {
    dataset: 'genomes',
    id: 'bed-genes',
    title: 'BED12 genes',
    summary:
      'Generated genes as BED12 on a 60 kb chromosome: a block an exon, each block a piece of its gene (one item, hovered and selected whole), the strand on each.',
    build: genomeBed.genes,
    registry: genomeBed.registry,
  },
  {
    dataset: 'genomes',
    id: 'bed-misread',
    title: 'A misread BED',
    summary:
      'A BED6 with faults: an end before its start, a position that is no number (both left out), a strand and a score that are not ones (read as missing), a line split on spaces and one with an extra column; each said once with its line number, the rest drawn.',
    build: genomeBed.misread,
    registry: genomeBed.registry,
  },
  {
    dataset: 'genomes',
    id: 'gff-tp53',
    title: 'TP53 from Ensembl’s GFF3',
    summary:
      'TP53’s region of Ensembl GRCh38 release 116’s GFF3 (the recorded slice, gzip), read whole by the opt-in `gff` runtime: its transcripts as gene models (exons as boxes, the CDS thick, the UTRs thin, a hairline with chevrons for the minus strand) and, from the same file with `features: genes`, one model a gene (TP53 and WRAP53’s start). chr17’s length from the file’s `##sequence-region`. Hover a transcript: every exon lights.',
    build: genomeGff.tp53Gff3,
    registry: genomeGff.registry,
  },
  {
    dataset: 'genomes',
    id: 'gtf-tp53',
    title: 'TP53 from GENCODE’s GTF',
    summary:
      'TP53’s region of GENCODE v50’s GTF (the recorded slice, gzip): the same gene models from GTF’s `transcript_id` and `gene_id`, the thick part to the stop codon GTF writes after the CDS; chr17’s length from an inline chrom.sizes.',
    build: genomeGff.tp53Gtf,
    registry: genomeGff.registry,
  },
  {
    dataset: 'genomes',
    id: 'gff-misread',
    title: 'A misread GFF3',
    summary:
      'A GFF3 with faults: an end before its start, a position that is no number and a short line (left out), an exon whose parent no line has, two transcripts each other’s parent (kept flat), a gene with a CDS and no transcript (read as one, said by type); each said once with its line number, the rest drawn as gene models.',
    build: genomeGff.misread,
    registry: genomeGff.registry,
  },
  {
    dataset: 'genomes',
    id: 'vcf-clinvar',
    title: 'TP53’s variants from ClinVar',
    summary:
      'TP53’s region of NCBI ClinVar’s GRCh38 sites (the recorded slice: BGZF, read by region at level 1 through the tabix index the source names) by the opt-in `vcf` runtime: a row an ALT allele, a mark at its REF extent coloured by its category (SNV, deletion, insertion, indel, MNV). ClinVar names chr17 `17`; its length from an inline chrom.sizes. Hover an allele: the card names it by ClinVar’s id with its track and extent; click selects its REF extent, a double-click focuses it.',
    build: genomeVcf.tp53Clinvar,
    registry: genomeVcf.registry,
  },
  {
    dataset: 'genomes',
    id: 'vcf-zoom-in',
    title: 'ClinVar’s chr17 at the fit',
    summary:
      'The same ClinVar file over the whole of chr17 (83 Mb): more level-1 tiles than a view reads, so none is asked and the track stays the band drawn while data is waiting (no text, no issue in the status line; no estimate drawn). Only TP53’s pages are recorded: zoom into TP53 (7.67–7.69 Mb) to read its variants.',
    build: genomeVcf.clinvarZoomIn,
    registry: genomeVcf.registry,
  },
  {
    dataset: 'genomes',
    id: 'vcf-long-deletion',
    title: 'A long deletion over SNVs',
    summary:
      'A generated VCF on a 20 kb chromosome (read whole; its length from its `##contig`): SNVs, multi-allelic sites (a row an allele), insertions, deletions and MNVs, and two long deletions, a symbolic `<DEL>` to its `END` and a 600-base sequence-resolved one, each one item drawn as a bar across the SNVs it covers.',
    build: genomeVcf.longDeletion,
    registry: genomeVcf.registry,
  },
  {
    dataset: 'genomes',
    id: 'vcf-misread',
    title: 'A misread VCF',
    summary:
      'A VCF with faults: a position that is no number, no REF and a short line (left out), a QUAL and a depth that are not numbers (read as missing), an undeclared INFO key, an `*` allele (no item), a record running past its contig (clipped) and one past it (left out), a header line without its Type; each said once with its line number, the rest drawn.',
    build: genomeVcf.misread,
    registry: genomeVcf.registry,
  },
];
