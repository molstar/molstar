/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// NCBI's genetic codes and the three-letter residue codes (seqspec:model#molecules, seqspec:model#codon-maps). Constant
// tables read without allocating.

/**
 * NCBI's genetic codes, copied from `gc.prt` Version 4.6 (lib/seqstar/core/_test/fixtures/gc.prt): per table id, its 64
 * amino acids (`ncbieaa`) then its start codons (`sncbieaa`, `M` a start), codons in TCAG order of their first,
 * second and third base. Twin: `GENETIC_CODES` in the builder's `helpers.ts`; a test keeps them equal.
 */
export const GC_PRT: Readonly<Record<number, string>> = {
  1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG---M------**--*----M---------------M----------------------------",
  2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG----------**--------------------MMMM----------**---M------------",
  3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------**----------------------MM---------------M------------",
  4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG--MM------**-------M------------MMMM---------------M------------",
  5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG---M------**--------------------MMMM---------------M------------",
  6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG--------------*--------------------M----------------------------",
  9: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG----------**-----------------------M---------------M------------",
  10: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------**-----------------------M----------------------------",
  11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG---M------**--*----M------------MMMM---------------M------------",
  12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------**--*----M---------------M----------------------------",
  13: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG---M------**----------------------MM---------------M------------",
  14: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG-----------*-----------------------M----------------------------",
  15: "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------*---*--------------------M----------------------------",
  16: "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------*---*--------------------M----------------------------",
  21: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG----------**-----------------------M---------------M------------",
  22: "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG------*---*---*--------------------M----------------------------",
  23: "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG--*-------**--*-----------------M--M---------------M------------",
  24: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG---M------**-------M---------------M---------------M------------",
  25: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG---M------**-----------------------M---------------M------------",
  26: "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------**--*----M---------------M----------------------------",
  27: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG--------------*--------------------M----------------------------",
  28: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------**--*--------------------M----------------------------",
  29: "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG--------------*--------------------M----------------------------",
  30: "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG--------------*--------------------M----------------------------",
  31: "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG----------**-----------------------M----------------------------",
  32: "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG---M------*---*----M------------MMMM---------------M------------",
  33: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG---M-------*-------M---------------M---------------M------------",
};

const X = 88;
/** A base's place in TCAG order by character code (T and U alike, either case); -1: an ambiguity code or other. */
const BASE = new Int8Array(128).fill(-1);
/** A base's complement by character code (T and U to A); others unchanged: an ambiguity code still translates to X. */
const COMPLEMENT = Uint8Array.from({ length: 128 }, (_, c) => c);
for (const [k, bases] of ["TUtu", "Cc", "Aa", "Gg"].entries()) for (const b of bases) BASE[b.charCodeAt(0)] = k;
for (const [b, c] of ["AT", "TA", "UA", "CG", "GC", "at", "ta", "ua", "cg", "gc"])
  COMPLEMENT[b.charCodeAt(0)] = c.charCodeAt(0);
const AAS: (Uint8Array | undefined)[] = [];
const STARTS: (Uint8Array | undefined)[] = [];
for (const [id, s] of Object.entries(GC_PRT)) {
  const codes = Uint8Array.from(s, (c) => c.charCodeAt(0));
  AAS[Number(id)] = codes.subarray(0, 64);
  STARTS[Number(id)] = codes.subarray(64);
}

/** A codon's index in TCAG order; -1 when a base is not one of T, U, C, A, G. */
function codonOf(a: number, b: number, c: number): number {
  const x = a < 128 ? BASE[a] : -1;
  const y = b < 128 ? BASE[b] : -1;
  const z = c < 128 ? BASE[c] : -1;
  return x < 0 || y < 0 || z < 0 ? -1 : (x << 4) | (y << 2) | z;
}

export const GeneticCode = {
  /** Whether NCBI defines the table. */
  has: (id: number): boolean => AAS[id] !== undefined,
  /** Three bases' character codes to a residue's (`*` a stop); an ambiguity code gives X; -1: no such table. */
  translate(id: number, a: number, b: number, c: number): number {
    const t = AAS[id];
    if (t === undefined) return -1;
    const k = codonOf(a, b, c);
    return k < 0 ? X : t[k];
  },
  /** Whether the codon is a start codon of the table. */
  isStart(id: number, a: number, b: number, c: number): boolean {
    const k = codonOf(a, b, c);
    return k >= 0 && STARTS[id]?.[k] === 77;
  },
  /** A base's complement (character codes): A↔T, C↔G, U → A; others unchanged. */
  complement: (c: number): number => (c < 128 ? COMPLEMENT[c] : c),
};

/** IUPAC's three-letter codes, and `Sec`, `Pyl`, `Asx`, `Glx`, `Xle`, `Xaa`, `Ter`. */
const THREE: Readonly<Record<string, string>> = Object.fromEntries(
  "A Ala R Arg N Asn D Asp C Cys Q Gln E Glu G Gly H His I Ile L Leu K Lys M Met F Phe P Pro S Ser T Thr W Trp Y Tyr V Val U Sec O Pyl B Asx Z Glx J Xle X Xaa * Ter"
    .split(" ")
    .flatMap((x, k, all) => (k % 2 ? [] : [[x, all[k + 1]]])),
);
/** A protein residue's three-letter code by its one-letter character code; `Xaa` for any other. For drawing (edges). */
export const threeLetter = (code: number): string => THREE[String.fromCharCode(code)] ?? "Xaa";
