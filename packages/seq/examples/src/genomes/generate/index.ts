/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MsgPackCodec } from '@molstar/core/repr/codecs/msgpack';
import { encodeSeq, type SeqData, seqToNode } from '@molstar/seqx-core/data';
import { generatedGenesBed, generatedGenotypesVcf } from '../files.js';
import {
  aliasTable,
  chromSizes,
  DefaultToyLengths,
  fasta,
  HlaName,
  type ToyAssembly,
  toyA,
  toyALetters,
} from './assemblies.js';
import { type AlignmentTruth, alignToyB, dotGrid, sidedPairs } from './alignment.js';
import { indexed } from './bgzf.js';
import { contactMap } from './contacts.js';
import { genesGff3, genesGtf, gffExtent, type ToyGene, toyGenes } from './genes.js';
import { type LiftoverTruth, liftoverBlocks, liftoverRuns, liftoverTruths, toyA2 } from './liftover.js';
import { cpgMethylation, geneMethylation } from './methylation.js';
import { bedExtent, bedGraph } from './signal.js';
import { GenomesSeed } from './streams.js';
import { SampleCount, SiteCount, type VariantsTruth, variantsVcf, vcfExtent } from './variants.js';

// The Genomes dataset's generator (examples.md, Genomes, The generator: outline): `generateGenomes(options)` returns
// every generated file of `fixtures/genomes/` as bytes (the record script writes them, a test compares them with the
// committed ones byte for byte) and the truths the views' tests read. Steps 1 to 6 and 8 here; step 7 (the ring's
// links) is the ring extension's. The VCF feature's `generated-genotypes.vcf` and the BED feature's
// `generated-genes.bed` are folded in from their generators, their bytes unchanged.

export interface GenomesOptions {
  /** `genomes.seed`: one fixed value by default. */
  readonly seed?: number;
  /** `genomes.toy.lengths`: toy A's chr1 to chr3 (300,000, 200,000, 120,000). */
  readonly lengths?: readonly [number, number, number];
  /** `genomes.samples`: the VCF's samples (40). */
  readonly samples?: number;
  /** `genomes.sites`: the VCF's sites (2,000). */
  readonly sites?: number;
}

/** The dataset's files, by their path under the package (`fixtures/genomes/…`). */
export const GenomesFiles = {
  sizesA: 'fixtures/genomes/toy-a.chrom.sizes',
  sizesB: 'fixtures/genomes/toy-b.chrom.sizes',
  sizesA2: 'fixtures/genomes/toy-a2.chrom.sizes',
  fasta: 'fixtures/genomes/toy-a.small.fa',
  aliases: 'fixtures/genomes/toy-a.aliases.seq',
  alignment: 'fixtures/genomes/toy-a-b.alignment.seq',
  dots: 'fixtures/genomes/toy-a-b.dots.seq.msgpack',
  liftover: 'fixtures/genomes/toy-a-a2.liftover.seq',
  gff3: 'fixtures/genomes/toy-a.genes.gff3.gz',
  gff3Index: 'fixtures/genomes/toy-a.genes.gff3.gz.tbi',
  gtf: 'fixtures/genomes/toy-a.genes.gtf.gz',
  gtfIndex: 'fixtures/genomes/toy-a.genes.gtf.gz.tbi',
  signal: 'fixtures/genomes/toy-a.signal.bedgraph.gz',
  signalIndex: 'fixtures/genomes/toy-a.signal.bedgraph.gz.tbi',
  vcf: 'fixtures/genomes/toy-a.variants.vcf.gz',
  vcfIndex: 'fixtures/genomes/toy-a.variants.vcf.gz.tbi',
  methylation: 'fixtures/genomes/toy-a.methylation.seq.msgpack',
  geneMethylation: 'fixtures/genomes/toy-a.gene-methylation.seq.msgpack',
  contacts: 'fixtures/genomes/toy-a.hic.seq.msgpack',
  genesBed: 'fixtures/genomes/generated-genes.bed',
  genotypes: 'fixtures/genomes/generated-genotypes.vcf',
} as const;

/** Each file's estimated size (examples.md, Genomes, Sizes), bytes: the generator's test asserts each under its own. */
export const EstimatedBytes: Readonly<Record<keyof typeof GenomesFiles, number>> = {
  sizesA: 1_000,
  sizesB: 1_000,
  sizesA2: 1_000,
  fasta: 22_000,
  aliases: 1_000,
  alignment: 25_000,
  dots: 60_000,
  liftover: 25_000,
  gff3: 60_000,
  gff3Index: 2_000,
  gtf: 60_000,
  gtfIndex: 2_000,
  signal: 120_000,
  signalIndex: 4_000,
  vcf: 60_000,
  vcfIndex: 2_000,
  methylation: 70_000,
  geneMethylation: 70_000,
  contacts: 60_000,
  genesBed: 2_000,
  genotypes: 25_000,
};

export interface GenomesTruth {
  readonly a: ToyAssembly;
  readonly a2: ToyAssembly;
  readonly alignment: AlignmentTruth;
  readonly genes: readonly ToyGene[];
  readonly liftover: readonly LiftoverTruth[];
  readonly variants: Omit<VariantsTruth, 'text'>;
  readonly cpgs: { readonly positions: readonly number[]; readonly missing: number };
  readonly lettersA: ReadonlyMap<string, string>;
}

export interface Genomes {
  /** Every generated file's bytes, by its path under the package. */
  readonly files: ReadonlyMap<string, Uint8Array>;
  readonly truth: GenomesTruth;
}

const enc = new TextEncoder();
/** A `seq` file in MessagePack (its typed columns as bytes: the numeric files, a third of their JSON). */
const msgpack = (data: SeqData) => MsgPackCodec.encode(seqToNode(data));
const names = (...as: ToyAssembly[]) => as.flatMap((a) => a.sequences.map((s) => s.name));

/** The dataset, generated: the same bytes on every run for the same options. */
export async function generateGenomes(options: GenomesOptions = {}): Promise<Genomes> {
  const seed = options.seed ?? GenomesSeed;
  const a = toyA(options.lengths ?? DefaultToyLengths);
  const lettersA = toyALetters(a, seed);
  const alignment = alignToyB(a, lettersA, seed);
  const b = alignment.b;
  const a2 = toyA2(a);
  const runs = liftoverRuns(a);
  const genes = toyGenes(a, seed);
  const variants = variantsVcf(a, lettersA, seed, options.samples ?? SampleCount, options.sites ?? SiteCount);
  const cpgs = cpgMethylation(a, lettersA, seed);
  const files = new Map<string, Uint8Array>();
  const put = (k: keyof typeof GenomesFiles, v: string | Uint8Array) =>
    files.set(GenomesFiles[k], typeof v === 'string' ? enc.encode(v) : v);

  put('sizesA', chromSizes(a));
  put('sizesB', chromSizes(b));
  put('sizesA2', chromSizes(a2));
  put('fasta', fasta(lettersA, ['chrM', HlaName]));
  put('aliases', encodeSeq(aliasTable(a)));
  put('alignment', encodeSeq(sidedPairs(a, b, alignment.blocks)));
  put('dots', msgpack(dotGrid(a, b, alignment.blocks)));
  put('liftover', encodeSeq(sidedPairs(a, a2, liftoverBlocks(runs))));
  const generic = { format: 0, zeroBased: false, meta: '#' } as const;
  const gff3 = await indexed(genesGff3(a, genes), { ...generic, columns: [1, 4, 5], names: names(a) }, gffExtent);
  put('gff3', gff3.file);
  put('gff3Index', gff3.index);
  const gtf = await indexed(genesGtf(a, genes), { ...generic, columns: [1, 4, 5], names: names(a) }, gffExtent);
  put('gtf', gtf.file);
  put('gtfIndex', gtf.index);
  const signal = await indexed(
    bedGraph([a], seed),
    { format: 0, zeroBased: true, meta: '#', columns: [1, 2, 3], names: names(a) },
    bedExtent,
  );
  put('signal', signal.file);
  put('signalIndex', signal.index);
  const vcf = await indexed(
    variants.text,
    { format: 2, zeroBased: false, meta: '#', columns: [1, 2, 0], names: names(a) },
    vcfExtent,
  );
  put('vcf', vcf.file);
  put('vcfIndex', vcf.index);
  put('methylation', msgpack(cpgs.data));
  put('geneMethylation', msgpack(geneMethylation(a, genes, seed)));
  put('contacts', msgpack(contactMap(a, seed)));
  put('genesBed', generatedGenesBed());
  put('genotypes', generatedGenotypesVcf());

  return {
    files,
    truth: {
      a,
      a2,
      alignment,
      genes,
      liftover: liftoverTruths(runs),
      variants: { sites: variants.sites },
      cpgs: { positions: cpgs.positions, missing: cpgs.missing },
      lettersA,
    },
  };
}
