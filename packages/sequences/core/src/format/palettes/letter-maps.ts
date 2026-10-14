/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The kernel's letter maps (seqspec:format#layer-types, `letter_map`): tables from a stored letter to the one shown,
// named data beside the letter color dictionaries. Letters not listed are kept. `complement` is the IUPAC nucleotide
// complement (ambiguity codes included); it is an instance of a general table, nothing in the core knows what a strand
// is: the entry lists the alphabets it is for, and resolution keeps any other alphabet's letters as stored.

import type { LetterMapEntry } from "../registry.ts";

/** Each two-letter word `xy` maps x to y, in both cases. */
const words = (s: string): Record<string, string> =>
  Object.fromEntries(
    s.split(" ").flatMap(([x, y]) => [
      [x, y],
      [x.toLowerCase(), y.toLowerCase()],
    ]),
  );
const IUPAC = "CG GC RY YR KM MK BV VB DH HD";

export const KernelLetterMaps: readonly LetterMapEntry[] = [
  { name: "complement", map: words(`AT TA UA ${IUPAC}`), alphabets: ["dna", "rna"] },
  { name: "complement_rna", map: words(`AU UA TA ${IUPAC}`), alphabets: ["dna", "rna"] },
  { name: "lower", map: Object.fromEntries([..."ABCDEFGHIJKLMNOPQRSTUVWXYZ"].map((c) => [c, c.toLowerCase()])) },
];

/** A character a letter map can name: printable ASCII (a space hides a letter). */
export const isMapLetter = (s: string) => s.length === 1 && s.charCodeAt(0) >= 32 && s.charCodeAt(0) < 127;

/**
 * A resolved `letter_map` as the 128-entry table a draw reads (`table[code]`: the code shown); null when it changes
 * nothing. Built once per prepared layer, so the per-letter lookup is one flat read.
 */
export function letterMapTable(map: Readonly<Record<string, string>> | null | undefined): Uint8Array | null {
  let table: Uint8Array | null = null;
  for (const [k, v] of Object.entries(map ?? {})) {
    if (!isMapLetter(k) || !isMapLetter(v) || k === v) continue;
    table ??= Uint8Array.from({ length: 128 }, (_, i) => i);
    table[k.charCodeAt(0)] = v.charCodeAt(0);
  }
  return table;
}
