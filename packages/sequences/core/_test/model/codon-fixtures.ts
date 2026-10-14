/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Case 14: a synthetic gene fragment and the protein its spliced CDS codes.

export const GENE =
  "GAATTCGCCACC" +
  "ATGGCTTCCAAGGAACTGACCCCACGTGACGGCC" +
  "GTAAGTCTTTCCTAACTTCTGTTTCCACAG" +
  "ATTGGTACTTCAACAGCCAGGTGGAGTAA" +
  "GCGGCCGCATAAAGT";
export const PROTEIN = "MASKELTPRDGHWYFNSQVE";
export const CDS: [number, number][] = [
  [12, 46],
  [76, 105],
];
