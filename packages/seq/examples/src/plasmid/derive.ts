/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Plasmid's derived data (examples.md, Plasmid): pUC19c's features found in its sequence, none taken from another
// annotation, its restriction sites computed on the cyclic sequence (a site across the origin found), the rotated
// variant, the reverse complement, and the GFF3, GenBank and `seq` texts of the committed fixtures. Plain functions of the
// letters with no imports (the record script runs them in Node as they are); a test compares each fixture with them.

/** The record: INSDC L09137.2, 2,686 bp, circular (its own title: Cloning vector pUC19c, complete sequence). */
export const PlasmidAccession = 'L09137.2';
/** The sequence's name in the fixtures (a `seq` system, a GFF3 `seqid`). */
export const PlasmidName = 'pUC19c';
export const PlasmidLength = 2686;

/** Where the rotated variant's origin is on the published one: inside bla, so its CDS crosses it. */
export const RotatedOrigin = 2000;

/** A feature: 0-based `start`, `length`, so its end may be past the sequence's length (it crosses the origin). */
export interface PlasmidFeature {
  readonly id: string;
  readonly label: string;
  /** The SO term (the GFF3 type). */
  readonly type: string;
  readonly strand: '+' | '-';
  readonly start: number;
  readonly length: number;
}

/** The complement of each IUPAC nucleotide code (R-Y, K-M, B-V, D-H pair up; S, W, N are their own). */
const Complement: Readonly<Record<string, string>> = {
  A: 'T',
  C: 'G',
  G: 'C',
  T: 'A',
  U: 'A',
  R: 'Y',
  Y: 'R',
  K: 'M',
  M: 'K',
  B: 'V',
  V: 'B',
  D: 'H',
  H: 'D',
  S: 'S',
  W: 'W',
  N: 'N',
};

export function reverseComplement(letters: string): string {
  let s = '';
  for (let i = letters.length - 1; i >= 0; i--) s += Complement[letters[i]] ?? 'N';
  return s;
}

/** The letters of a construct whose origin is at `origin` of the published one. */
export const rotateLetters = (letters: string, origin: number) => letters.slice(origin) + letters.slice(0, origin);

/** A position of the published construct on one whose origin is at `origin`. */
export const rotate = (p: number, origin: number, length: number = PlasmidLength) =>
  (((p - origin) % length) + length) % length;

/** The standard genetic code (TCAG order). */
const Code = (() => {
  const aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  const out: Record<string, string> = {};
  let i = 0;
  for (const x of 'TCAG') for (const y of 'TCAG') for (const z of 'TCAG') out[x + y + z] = aa[i++];
  return out;
})();

export function translate(dna: string): string {
  let p = '';
  for (let i = 0; i + 2 < dna.length; i += 3) p += Code[dna.slice(i, i + 3)] ?? 'X'; // an ambiguous codon is X
  return p;
}

/** The text of `[start, start + length)` on a cyclic sequence. */
function cyclicSlice(letters: string, start: number, length: number): string {
  const end = start + length;
  return end <= letters.length
    ? letters.slice(start, end)
    : letters.slice(start) + letters.slice(0, end - letters.length);
}

/**
 * The open reading frame, on either strand of the cyclic sequence, that starts at an ATG and encodes a protein starting
 * with `prefix` (read to its first stop, the stop included in the feature). `start` is on the published strand.
 */
export function findOrf(
  letters: string,
  prefix: string,
): { start: number; length: number; strand: '+' | '-' } | undefined {
  const n = letters.length;
  for (const strand of ['+', '-'] as const) {
    const s = strand === '+' ? letters : reverseComplement(letters);
    const doubled = s + s;
    for (let at = 0; at < n; at++) {
      if (!doubled.startsWith('ATG', at)) continue;
      // the first residues first: most starts are not the protein's
      if (translate(doubled.slice(at, at + 3 * prefix.length)) !== prefix) continue;
      let k = at;
      let p = '';
      while (k < at + n) {
        const aa = Code[doubled.slice(k, k + 3)] ?? 'X';
        k += 3;
        p += aa;
        if (aa === '*') break;
      }
      if (!p.endsWith('*') || !p.startsWith(prefix)) continue;
      const length = k - at;
      // on the reverse complement a position `at` is `n - at - length` as a start on the published strand
      return { start: strand === '+' ? at : (((n - at - length) % n) + n) % n, length, strand };
    }
  }
  return undefined;
}

/** The first span, on either strand, from `first` to `last` (the signature's two ends), found in the cyclic sequence. */
export function findSpan(
  letters: string,
  first: string,
  last: string,
): { start: number; length: number; strand: '+' | '-' } | undefined {
  const n = letters.length;
  for (const strand of ['+', '-'] as const) {
    const s = strand === '+' ? letters : reverseComplement(letters);
    const doubled = s + s;
    const a = doubled.indexOf(first);
    if (a < 0 || a >= n) continue;
    const b = doubled.indexOf(last, a);
    if (b < 0) continue;
    const length = b + last.length - a;
    return { start: strand === '+' ? a : (((n - a - length) % n) + n) % n, length, strand };
  }
  return undefined;
}

/**
 * The pMB1 origin of replication (the ColE1-type origin of pBR322, Sutcliffe 1979, in pUC vectors, Yanisch-Perron et
 * al. 1985), located by the two ends of the stretch plasmid maps annotate as `ori` (589 bp in pUC19). Written here, the
 * ends found in the record on its minus strand and checked by the anchors.
 */
export const OriSignature = {
  first: 'TTGAGATCCTTTTTTTCTGCGCGTAATCTGC',
  last: 'GGGGGGCGGAGCCTATGGAAA',
} as const;

/** The multiple cloning site runs between the sites of its first enzyme and its last (the shorter way round). */
export const McsEnds = { first: 'EcoRI', last: 'HindIII' } as const;

/** The proteins whose open reading frames are the CDS features, by their first residues (TEM-1 β-lactamase; lacZα). */
export const BlaPrefix = 'MSIQHFRVAL';
export const LacZaPrefix = 'MTMITPSLHA';

/** The restriction enzymes of the sites layer: the cloning site's, then unique cutters elsewhere. All palindromic. */
export const Enzymes: readonly {
  readonly name: string;
  readonly site: string;
  /** Where the top strand is cut: bases of the site before the cut (EcoRI G^AATTC: 1). */
  readonly cut: number;
  readonly mcs: boolean;
}[] = [
  { name: 'EcoRI', site: 'GAATTC', cut: 1, mcs: true },
  { name: 'SacI', site: 'GAGCTC', cut: 5, mcs: true },
  { name: 'KpnI', site: 'GGTACC', cut: 5, mcs: true },
  { name: 'SmaI', site: 'CCCGGG', cut: 3, mcs: true },
  { name: 'BamHI', site: 'GGATCC', cut: 1, mcs: true },
  { name: 'XbaI', site: 'TCTAGA', cut: 1, mcs: true },
  { name: 'SalI', site: 'GTCGAC', cut: 1, mcs: true },
  { name: 'PstI', site: 'CTGCAG', cut: 5, mcs: true },
  { name: 'SphI', site: 'GCATGC', cut: 5, mcs: true },
  { name: 'HindIII', site: 'AAGCTT', cut: 1, mcs: true },
  { name: 'NdeI', site: 'CATATG', cut: 2, mcs: false },
  { name: 'ScaI', site: 'AGTACT', cut: 3, mcs: false },
  { name: 'AatII', site: 'GACGTC', cut: 5, mcs: false },
];

/** A recognition site: the enzyme, where it starts (0-based, on the published strand; its end may pass the length). */
export interface PlasmidSite {
  readonly enzyme: string;
  readonly start: number;
  readonly length: number;
  /** The position the top strand is cut before (0-based on the sequence: past the origin, wrapped). */
  readonly cut: number;
}

/** Every site of every enzyme on the cyclic sequence, a site across the origin included (all enzymes palindromic). */
export function findSites(letters: string): PlasmidSite[] {
  const n = letters.length;
  const out: PlasmidSite[] = [];
  for (const e of Enzymes) {
    const doubled = letters + letters.slice(0, e.site.length - 1);
    for (let i = doubled.indexOf(e.site); i >= 0 && i < n; i = doubled.indexOf(e.site, i + 1))
      out.push({ enzyme: e.name, start: i, length: e.site.length, cut: (i + e.cut) % n });
  }
  return out.sort((a, b) => a.start - b.start);
}

/**
 * The features found in the sequence (on the published origin): bla and lacZα as the open reading frames that encode
 * their proteins, the multiple cloning site from its first site (EcoRI) to its last (HindIII), the shorter way round, the origin by its
 * signature. Sorted by start.
 */
export function deriveFeatures(letters: string): PlasmidFeature[] {
  const out: PlasmidFeature[] = [];
  const bla = findOrf(letters, BlaPrefix);
  if (bla) out.push({ id: 'bla', label: 'bla (AmpR)', type: 'CDS', ...bla });
  const lacZa = findOrf(letters, LacZaPrefix);
  if (lacZa) out.push({ id: 'lacZa', label: 'lacZα', type: 'CDS', ...lacZa });
  const all = findSites(letters);
  const first = all.find((x) => x.enzyme === McsEnds.first);
  const last = all.find((x) => x.enzyme === McsEnds.last);
  if (first && last) {
    const n = letters.length;
    const mod = (x: number) => ((x % n) + n) % n;
    // the span between the two sites that holds them both: the shorter way round, whichever strand the record is on
    const forward = mod(last.start + last.length - first.start) || n;
    const back = mod(first.start + first.length - last.start) || n;
    out.push({
      id: 'mcs',
      label: 'MCS',
      type: 'misc_feature',
      strand: '+',
      start: forward <= back ? first.start : last.start,
      length: Math.min(forward, back),
    });
  }
  const ori = findSpan(letters, OriSignature.first, OriSignature.last);
  if (ori) out.push({ id: 'ori', label: 'ori', type: 'rep_origin', ...ori });
  return out.sort((a, b) => a.start - b.start);
}

/** The rotated variant's primer-binding site: generated, across its origin on the plus strand (EC-11). */
export function primerSite(origin: number = RotatedOrigin): PlasmidFeature {
  return {
    id: 'primer',
    label: 'primer site',
    type: 'primer_bind',
    strand: '+',
    start: origin - 30,
    length: 55,
  };
}

/** A feature on a construct whose origin is at `origin` of the published one (its length kept). */
export const rotateFeature = <T extends { start: number }>(f: T, origin: number, length = PlasmidLength): T => ({
  ...f,
  start: rotate(f.start, origin, length),
});

/** The text of a feature's whole letters, on its own strand. */
export function featureLetters(letters: string, f: { start: number; length: number; strand: '+' | '-' }): string {
  const s = cyclicSlice(letters, f.start, f.length);
  return f.strand === '+' ? s : reverseComplement(s);
}

const gffEscape = (v: string) =>
  v.replace(/[;=&,%\t\n]/g, (c) => `%${c.charCodeAt(0).toString(16).toUpperCase().padStart(2, '0')}`);

/** A GenBank location of `[start, start + length)` on a cyclic sequence: a join over the origin, `complement` on minus. */
function genbankLocation(f: { start: number; length: number; strand: '+' | '-' }, length: number): string {
  const end = f.start + f.length;
  const plain = end <= length ? `${f.start + 1}..${end}` : `join(${f.start + 1}..${length},1..${end - length})`;
  return f.strand === '-' ? `complement(${plain})` : plain;
}

/** A qualifier's lines (from column 22, to 79), a value cut where it reaches the edge. */
function qualifierLines(name: string, value: string): string[] {
  const indent = ' '.repeat(21);
  const text = `/${name}="${value.replaceAll('"', '""')}"`;
  const out: string[] = [];
  for (let k = 0; k < text.length; k += 58) out.push(indent + text.slice(k, k + 58));
  return out;
}

/**
 * The rotated variant as a GenBank flat file, as NCBI writes a record (LOCUS, DEFINITION, SOURCE, a REFERENCE and a
 * COMMENT, the feature table, ORIGIN's lowercase letters): `source`, then the features in the GFF3's order, a feature
 * across the origin a `join` over it (bla's on the minus strand `complement(join(…))`), each labelled as the GFF3's
 * `Name` (`/label`), a CDS's protein in `/translation`. The sites are no features (INSDC has no key for them).
 */
export function genbankText(
  letters: string,
  features: readonly PlasmidFeature[],
  length: number = PlasmidLength,
): string {
  const pad = (s: string, n: number) => s + ' '.repeat(Math.max(1, n - s.length));
  const lines = [
    `LOCUS       ${pad(PlasmidName, 16)}${String(length).padStart(11)} bp    DNA     circular SYN 19-NOV-2026`,
    `DEFINITION  Cloning vector ${PlasmidName}, turned to an origin at ${RotatedOrigin + 1}, inside bla; its features`,
    '            derived in code (generated).',
    'KEYWORDS    .',
    `SOURCE      Cloning vector ${PlasmidName}`,
    `  ORGANISM  Cloning vector ${PlasmidName}`,
    '            other sequences; artificial sequences; vectors.',
    `REFERENCE   1  (bases 1 to ${length})`,
    '  AUTHORS   Yanisch-Perron,C., Vieira,J. and Messing,J.',
    '  TITLE     Improved M13 phage cloning vectors and host strains: nucleotide',
    '            sequences of the M13mp18 and pUC19 vectors',
    '  JOURNAL   Gene 33 (1), 103-119 (1985)',
    `COMMENT     Generated by the Seq* examples from the letters of INSDC ${PlasmidAccession}.`,
    'FEATURES             Location/Qualifiers',
    `     source          1..${length}`,
    ...qualifierLines('organism', `Cloning vector ${PlasmidName}`),
    ...qualifierLines('mol_type', 'other DNA'),
  ];
  for (const f of features) {
    lines.push(`     ${pad(f.type, 16)}${genbankLocation(f, length)}`);
    if (f.type === 'CDS') lines.push(...qualifierLines('gene', f.id));
    lines.push(...qualifierLines('label', f.label));
    if (f.id === 'primer') lines.push(...qualifierLines('note', 'generated, across the origin'));
    if (f.type === 'CDS') {
      lines.push('                     /codon_start=1');
      lines.push(...qualifierLines('translation', translate(featureLetters(letters, f)).replace(/\*$/, '')));
    }
  }
  lines.push('ORIGIN      ');
  const lower = letters.toLowerCase();
  for (let at = 0; at < lower.length; at += 60) {
    const blocks: string[] = [];
    for (let k = at; k < Math.min(at + 60, lower.length); k += 10) blocks.push(lower.slice(k, k + 10));
    lines.push(`${String(at + 1).padStart(9)} ${blocks.join(' ')}`);
  }
  lines.push('//');
  return `${lines.join('\n')}\n`;
}

/**
 * The features and the sites as GFF3, as NCBI writes a ring's: a `region` row (`Is_circular`), a feature across the
 * origin one row ending past the length (the records runtime's cyclic wrap). Positions from 1, the last one named.
 */
export function gff3Text(
  features: readonly PlasmidFeature[],
  sites: readonly PlasmidSite[],
  length: number = PlasmidLength,
): string {
  const lines = ['##gff-version 3', `##sequence-region ${PlasmidName} 1 ${length}`];
  const row = (type: string, start: number, len: number, strand: string, attributes: string) =>
    [PlasmidName, 'seqx', type, start + 1, start + len, '.', strand, '.', attributes].join('\t');
  lines.push(row('region', 0, length, '+', `ID=region-${PlasmidName};Name=${PlasmidName};Is_circular=true`));
  for (const f of features) {
    const note = f.id === 'primer' ? ';note=generated, across the origin' : '';
    lines.push(row(f.type, f.start, f.length, f.strand, `ID=${f.id};Name=${gffEscape(f.label)}${note}`));
  }
  for (const s of sites)
    lines.push(row('restriction_enzyme_recognition_site', s.start, s.length, '+', `ID=${s.enzyme};Name=${s.enzyme}`));
  return `${lines.join('\n')}\n`;
}

/**
 * The `seq` file of the cyclic sequence: one system, one region over all of it at level 1 (what `encodeSeq` of
 * `inlineSeq` writes; a test compares the two), compact.
 */
export function seqText(name: string, letters: string): string {
  const system = { name, length: letters.length, alphabet: 'dna', cyclic: true };
  const shapes = ['letters'];
  const data = {
    kind: 'seq',
    params: {
      coverage: { systems: [system], shapes: [{ system: name, shape: 'letters' }], levels: [1] },
      pairs: [
        {
          request: { system: name, start: 0, end: letters.length, level: 1, shapes },
          answer: { ok: true, tables: { letters: { form: 'letters', start: 0, letters } } },
        },
      ],
    },
  };
  return `${JSON.stringify(data)}\n`;
}

/**
 * The fixtures' locations: the record as a `seq` file, the rotated variant's letters as one (the record turned to its
 * origin), and the rotated variant's features served as GFF3, and as GenBank with its letters.
 */
export const PlasmidFiles = {
  published: 'fixtures/plasmid/L09137.2.seq',
  rotated: 'fixtures/plasmid/pUC19c-rotated.seq',
  rotatedFeatures: 'fixtures/plasmid/pUC19c-rotated-features.gff3',
  rotatedGenbank: 'fixtures/plasmid/pUC19c-rotated.gb',
} as const;

/** The module the record script writes from the record: its letters, what the annotations are derived from in code. */
export function recordModule(letters: string): string {
  return `/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// INSDC ${PlasmidAccession}, Cloning vector ${PlasmidName}, complete sequence (${letters.length} bp, circular): the letters of
// \`fixtures/plasmid/${PlasmidAccession}.seq\`, written by \`scripts/record.ts\` (a test compares the two): what the
// features and sites are derived from in code, at once when a view is built. The views read their letters from the
// \`seq\` files.

export const PublishedLetters =
  '${letters}';
`;
}

/** The features and sites of the rotated variant, from the published letters: derived, turned, the primer site added. */
export function rotatedAnnotations(letters: string, origin: number = RotatedOrigin) {
  const turned = rotateLetters(letters, origin);
  const features = [
    ...deriveFeatures(letters).map((f) => rotateFeature(f, origin)),
    rotateFeature(primerSite(origin), origin),
  ];
  return { turned, features, sites: findSites(turned) };
}
