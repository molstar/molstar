/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The builder's own runtime values (seqspec:format#builders: no runtime dependencies). Each is a small copy; its twin
// is named beside it, and a change to either is made to both.

/** Twin: `SvsFormat` and `SvsVersion` in the format's `doc.ts`. */
export const SvsFormat = 'seq_view_spec';
export const SvsVersion = '1.0.0';

/** Twin: `NAMESPACED` and `isNamespaced` in the format's `schema/field-schema.ts`. */
const NAMESPACED = /^[a-z][a-z0-9_-]*(\.[a-z0-9_-]+)*:[A-Za-z0-9][A-Za-z0-9_.-]*$/;
export const isNamespaced = (name: unknown): name is string => typeof name === 'string' && NAMESPACED.test(name);

/** A CIGAR from two gapped rows of equal length (`-` a gap; a column gapped in both skipped). Twin: `cigarFromRows`
 *  in algorithms' `src/aligner.ts`. */
export function cigarFromRows(a: string, b: string): string {
  if (a.length !== b.length) throw new RangeError(`rows of different lengths: ${a.length} and ${b.length}`);
  let cigar = '';
  let op = '';
  let run = 0;
  for (let i = 0; i < a.length; i++) {
    const inA = a[i] !== '-';
    const inB = b[i] !== '-';
    if (!inA && !inB) continue;
    const o = inA && inB ? 'M' : inA ? 'I' : 'D';
    if (o !== op && run > 0) cigar += `${run}${op}`;
    run = o === op ? run + 1 : 1;
    op = o;
  }
  return run > 0 ? `${cigar}${run}${op}` : cigar;
}

/** Twin: `float32Number` in model's `src/common.ts` (why it is shortest there); a test keeps them equal. */
export function f32(x: number): number {
  if (x === 0 || !Number.isFinite(x)) return x + 0;
  const y = Number(x.toPrecision(7));
  if (Math.fround(y) === x) {
    const d = 0.999999 * 10 ** (Math.floor(Math.log10(Math.abs(y))) - 6);
    if (Math.fround(y + d) !== x && Math.fround(y - d) !== x) return y;
    for (let p = 1; p < 7; p++) if (Math.fround(Number(x.toPrecision(p))) === x) return Number(x.toPrecision(p));
    return y;
  }
  const z = Number(x.toPrecision(8));
  if (Math.fround(z) === x) return z;
  const w = Number(x.toPrecision(9));
  const u = Math.abs(x) === 2 ** Math.floor(Math.log2(Math.abs(x))) ? Number((2 * w - z).toPrecision(8)) : w;
  return Math.fround(u) === x ? u : w;
}

/** WUSS pairs: per column its partner, -1 unpaired; `()`, `<>`, `[]`, `{}` and `Aa` letter pairs pair. Unbalanced
 *  input throws. Twin: `wussPairs` (`Stockholm.pairs`) in model's `src/formats.ts`, which reports an issue instead. */
export function wussPairs(ss: string): Int32Array {
  const partner = new Int32Array(ss.length).fill(-1);
  const stacks: number[][] = Array.from({ length: 30 }, () => []);
  for (let i = 0; i < ss.length; i++) {
    const code = ss.charCodeAt(i);
    let open = '([{<'.indexOf(ss[i]);
    let close = ')]}>'.indexOf(ss[i]);
    if (code >= 65 && code <= 90) open = 4 + code - 65;
    else if (code >= 97 && code <= 122) close = 4 + code - 97;
    if (open >= 0) stacks[open].push(i);
    else if (close >= 0) {
      const j = stacks[close].pop();
      if (j === undefined) throw new Error(`basePairs: WUSS: "${ss[i]}" at column ${i} closes nothing.`);
      partner[i] = j;
      partner[j] = i;
    }
  }
  const left = stacks.findIndex((s) => s.length > 0);
  if (left >= 0) throw new Error(`basePairs: WUSS: the pair opened at column ${stacks[left][0]} is not closed.`);
  return partner;
}

/**
 * NCBI's genetic codes (`gc.prt` Version 4.6): per table id, its 64 amino acids then its starts (`M`), codons in TCAG
 * order. Twin: `GC_PRT` in the model's `genetic-code.ts`; a test in the format keeps them equal.
 */
export const GENETIC_CODES: Readonly<Record<number, string>> = {
  1: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG---M------**--*----M---------------M----------------------------',
  2: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG----------**--------------------MMMM----------**---M------------',
  3: 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------**----------------------MM---------------M------------',
  4: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG--MM------**-------M------------MMMM---------------M------------',
  5: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG---M------**--------------------MMMM---------------M------------',
  6: 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG--------------*--------------------M----------------------------',
  9: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG----------**-----------------------M---------------M------------',
  10: 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------**-----------------------M----------------------------',
  11: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG---M------**--*----M------------MMMM---------------M------------',
  12: 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------**--*----M---------------M----------------------------',
  13: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG---M------**----------------------MM---------------M------------',
  14: 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG-----------*-----------------------M----------------------------',
  15: 'FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------*---*--------------------M----------------------------',
  16: 'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------*---*--------------------M----------------------------',
  21: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG----------**-----------------------M---------------M------------',
  22: 'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG------*---*---*--------------------M----------------------------',
  23: 'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG--*-------**--*-----------------M--M---------------M------------',
  24: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG---M------**-------M---------------M---------------M------------',
  25: 'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG---M------**-----------------------M---------------M------------',
  26: 'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------**--*----M---------------M----------------------------',
  27: 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG--------------*--------------------M----------------------------',
  28: 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------**--*--------------------M----------------------------',
  29: 'FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG--------------*--------------------M----------------------------',
  30: 'FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG--------------*--------------------M----------------------------',
  31: 'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------**-----------------------M----------------------------',
  32: 'FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG---M------*---*----M------------MMMM---------------M------------',
  33: 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG---M-------*-------M---------------M---------------M------------',
};

const COMPLEMENT: Readonly<Record<string, string>> = { A: 'T', T: 'A', U: 'A', C: 'G', G: 'C' };
const AMBIGUITY = 'RYSWKMBDHV';
const upper = (c: string) => (c >= 'a' && c <= 'z' ? c.toUpperCase() : c);

/**
 * Why `dna` is not DNA or RNA as the model reads one, or null. Twin: `LETTERS` and `detectAlphabet` in model's
 * `src/molecule.ts` (letters by ASCII upper case; given `alphabet`: every letter is one of it; else nucleotide when
 * every letter is a nucleotide code and at least 90% are `ACGTUN`, RNA when `U` is there and `T` is not).
 */
export function notNucleotide(dna: string, alphabet: unknown): string | null {
  const n: Record<string, number> = { T: 0, U: 0 };
  let plain = 0;
  let other = 0;
  for (const ch of dna) {
    const c = upper(ch);
    if ('ACGTUN'.includes(c)) plain++;
    else if (!AMBIGUITY.includes(c)) other++;
    if (c in n) n[c]++;
  }
  if (alphabet !== undefined && alphabet !== null && alphabet !== 'dna' && alphabet !== 'rna')
    return `its alphabet is ${String(alphabet)}`;
  const kind = alphabet ?? (!other && plain * 10 >= dna.length * 9 && dna ? (n.U && !n.T ? 'rna' : 'dna') : 'protein');
  if (kind === 'protein') return 'its letters are not nucleotide codes';
  if (other || n[kind === 'dna' ? 'U' : 'T']) return `a letter is not one of ${kind}'s`;
  return null;
}
/** A codon's residue (`*` a stop; X when a base is an ambiguity code) and whether it is a start codon. */
export function codonOf(table: string, codon: string): { residue: string; start: boolean } {
  const k = [...codon].reduce((x, b) => {
    const i = 'TCAG'.indexOf(b === 'U' ? 'T' : b);
    return x < 0 || i < 0 ? -1 : x * 4 + i;
  }, 0);
  return k < 0 ? { residue: 'X', start: false } : { residue: table[k], start: table[64 + k] === 'M' };
}

/**
 * A coding range translated: the first codon `M` when a start of the table, a final stop left out, an inner one `*`;
 * with `cds` null, from `start` to the first stop codon or the last whole codon. Throws on bad input.
 */
export function translateDna(dna: string, cds: [number, number][] | null, start: number, strand: string, code: number) {
  const table = GENETIC_CODES[code];
  if (!table) throw new RangeError(`translate: NCBI defines no genetic code ${code}.`);
  if (strand === '-' && !cds)
    throw new RangeError('translate: on strand - give `cds`: its ranges fix where reading starts.');
  const s = dna.toUpperCase();
  cds?.forEach(([a, b], k) => {
    // The codon mapping's order rule (twin: `codons` in model's `src/mapping.ts`).
    if (k > 0 && (strand === '-' ? b > cds[k - 1][0] : a < cds[k - 1][1])) {
      const order = strand === '-' ? 'descending (strand -)' : 'ascending';
      throw new RangeError(`translate: range ${k} overlaps range ${k - 1} or is out of ${order} order.`);
    }
  });
  const read = (a: number, b: number) => {
    if (!(Number.isInteger(a) && Number.isInteger(b) && a >= 0 && a < b && b <= s.length))
      throw new RangeError(`translate: [${a}, ${b}) is outside the sequence.`);
    return strand === '-'
      ? [...s.slice(a, b)]
          .reverse()
          .map((c) => COMPLEMENT[c] ?? c)
          .join('')
      : s.slice(a, b);
  };
  const nts = cds ? cds.map(([a, b]) => read(a, b)).join('') : read(start, s.length);
  if (cds && nts.length % 3)
    throw new RangeError(`translate: the CDS has ${nts.length} nucleotides, not a multiple of 3.`);
  let out = '';
  for (let i = 0; i + 3 <= nts.length; i += 3) {
    const { residue, start: first } = codonOf(table, nts.slice(i, i + 3));
    if (residue === '*' && (!cds || i + 3 === nts.length)) break;
    out += i === 0 && first ? 'M' : residue;
  }
  if (!out) throw new RangeError('translate: the translation is empty: its first codon is a stop.');
  return out;
}
