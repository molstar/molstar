/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { encodeSeq, type SeqData, inlineSeq } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';

// The p53 dataset: human TP53 (UniProt P04637), 393 residues, its domains and a conservation signal per residue, and a
// toy gene it is placed on (minus strand, three coding exons). The protein is the real one; the domains' bounds are
// illustrative; the conservation and the gene's introns and flanks are generated from a seed and never stored; the
// gene's coding exons encode the protein (one codon a residue, the same for each amino acid) but for one known variant,
// R175H, so a view checking codons shows the conflict where it is. The gene's letters are also committed as a small
// file (`fixtures/p53/tp53-gene.seq`) that the views read through the source protocol, the way data that is not inline
// is read.

/** Human p53 (UniProt P04637), 393 residues. */
export const P53 =
  'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD';

/** A seeded generator (a 32-bit LCG): generated data is made from a seed, never stored. */
export function seeded(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 0x100000000;
  };
}

/** Conservation: a value a residue in [0, 1], two missing (NaN). */
export function conservation(): Float32Array {
  const r = seeded(53);
  const v = Float32Array.from({ length: P53.length }, () => Math.round(r() * 100) / 100);
  v[100] = Number.NaN;
  v[101] = Number.NaN;
  return v;
}

/** The domains: TAD, proline-rich, DNA-binding, nuclear localization, tetramerization, C-terminal. */
export const Domains = {
  start: new Int32Array([0, 63, 93, 318, 323, 363]),
  end: new Int32Array([61, 94, 292, 360, 355, 393]),
  id: ['TAD', 'PRR', 'DBD', 'NLS', 'TET', 'CTD'],
  label: ['Transactivation', 'Proline-rich', 'DNA-binding', 'Nuclear localization', 'Tetramerization', 'C-terminal'],
  category: ['activation', 'linker', 'binding', 'signal', 'oligomer', 'regulatory'],
};

/**
 * p53's secondary structure: helices and strands of the core domain and the tetramerization helix (bounds illustrative,
 * as the domains'), the `category` the class.
 */
export const Secondary = {
  start: new Int32Array([16, 100, 111, 124, 135, 146, 156, 171, 203, 213, 227, 233, 251, 264, 277, 334]),
  end: new Int32Array([25, 110, 120, 133, 141, 152, 166, 178, 210, 224, 232, 240, 260, 273, 287, 343]),
  id: ['H1', 'S1', 'S2', 'S3', 'H2', 'S4', 'S5', 'H3', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'H4', 'H5'],
  category: [
    'helix',
    'strand',
    'strand',
    'strand',
    'helix',
    'strand',
    'strand',
    'helix',
    'strand',
    'strand',
    'strand',
    'strand',
    'strand',
    'strand',
    'helix',
    'helix',
  ],
};

/** p53's well-known variants (position, change, consequence), two or three at the hot spots R248 and R273. */
export const Variants = {
  start: new Int32Array([174, 195, 212, 219, 244, 247, 247, 248, 272, 272, 281, 305]),
  end: new Int32Array([175, 196, 213, 220, 245, 248, 248, 249, 273, 273, 282, 306]),
  id: ['R175H', 'R196*', 'R213*', 'Y220C', 'G245S', 'R248Q', 'R248W', 'R249S', 'R273H', 'R273C', 'R282W', 'R306*'],
  label: ['R175H', 'R196*', 'R213*', 'Y220C', 'G245S', 'R248Q', 'R248W', 'R249S', 'R273H', 'R273C', 'R282W', 'R306*'],
  category: [
    'missense',
    'nonsense',
    'nonsense',
    'missense',
    'missense',
    'missense',
    'missense',
    'missense',
    'missense',
    'missense',
    'missense',
    'nonsense',
  ],
};

type SpecBuilder = ReturnType<typeof createSeqSpec>;

/** The secondary structure as an annotation of the protein. */
export function addSecondary(b: SpecBuilder, protein: Parameters<SpecBuilder['data']['ranges']>[1]['target']) {
  return b.data.ranges('secondary', { target: protein, rows: Secondary });
}

/** The variants as an annotation of the protein. */
export function addVariants(b: SpecBuilder, protein: Parameters<SpecBuilder['data']['ranges']>[1]['target']) {
  return b.data.ranges('variants', { target: protein, rows: Variants });
}

/** The domains as an annotation of the protein (what a view places on a track). */
export function addDomains(b: SpecBuilder, protein: Parameters<SpecBuilder['data']['ranges']>[1]['target']) {
  return b.data.ranges('domains', { target: protein, rows: Domains });
}

/** The toy gene's exons (minus strand, genomic order): 300 + 420 + 462 = 1182 nucleotides, 393 codons and a stop. */
export const GeneExons = { start: [100, 480, 960], end: [400, 900, 1422] };

export const GeneLength = 1500;

/** A codon for each amino acid (the one human genes use most), and the stop. */
const Codons: Readonly<Record<string, string>> = {
  A: 'GCC',
  R: 'CGG',
  N: 'AAC',
  D: 'GAC',
  C: 'TGC',
  Q: 'CAG',
  E: 'GAG',
  G: 'GGC',
  H: 'CAC',
  I: 'ATC',
  L: 'CTG',
  K: 'AAG',
  M: 'ATG',
  F: 'TTC',
  P: 'CCC',
  S: 'AGC',
  T: 'ACC',
  W: 'TGG',
  Y: 'TAC',
  V: 'GTG',
  '*': 'TGA',
};

/** The variant the toy gene carries: R175H (residue 174 from 0, a hotspot), its codon His where the protein says Arg. */
export const GeneVariant = { residue: 174, from: 'R', to: 'H' } as const;

/** The coding sequence, 5' to 3': p53 back-translated, the variant's codon, then the stop (1182 nucleotides). */
export function geneCoding(): string {
  let s = '';
  for (let i = 0; i < P53.length; i++) s += Codons[i === GeneVariant.residue ? GeneVariant.to : P53[i]];
  return `${s}${Codons['*']}`;
}

const Complement: Readonly<Record<string, string>> = { A: 'T', C: 'G', G: 'C', T: 'A' };

/**
 * The toy gene's letters: the flanks and introns from a seed; the coding sequence on the minus strand, read from the
 * last exon's end down, complemented.
 */
export function geneLetters(): string {
  const r = seeded(17);
  const s: string[] = [];
  for (let i = 0; i < GeneLength; i++) s.push('ACGT'[Math.floor(r() * 4)]);
  const coding = geneCoding();
  let c = 0;
  for (let k = GeneExons.start.length - 1; k >= 0; k--)
    for (let g = GeneExons.end[k] - 1; g >= GeneExons.start[k]; g--) s[g] = Complement[coding[c++]];
  return s.join('');
}

/** The committed file of the gene's letters: the source a view names by location. */
export const GeneFile = 'fixtures/p53/tp53-gene.seq';

/** The gene's letters as the committed file's text (a source of one region over the whole system, level 1). */
export function geneSeqText(): string {
  const data: SeqData = inlineSeq(
    { name: 'tp53', length: GeneLength, alphabet: 'dna' },
    { letters: { form: 'letters', start: 0, letters: geneLetters() } },
  );
  return `${encodeSeq(data)}\n`;
}
