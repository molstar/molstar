/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * The `clustal_x` theme: Clustal X coloring by a column's composition, as Jalview's `ClustalxColourScheme`
 * (seqspec:format#themes).
 *
 * Written from Jalview's documented rules (https://www.jalview.org/help/html/colourSchemes/clustal.html, the table of
 * residue groups and thresholds; Jalview is GPL-3.0, only its published rules are used here, no code). Each residue
 * is colored by the first of its rules whose criterion holds in the column: a share of a residue group above a
 * threshold, or (a list of letters) the share of any one of them above it. Shares are of the column's letters; gaps
 * (anything not a letter) are left out. Lower case counts as upper case. The colors are Clustal's groups as MVS has
 * them (`AminoGroupColors`). Protein only (a viewer uses `letter` on nucleotides).
 *
 * Choices where the documented table is open: thresholds are strict (`> 60%`); the cysteine rule (`> 85%` C, pink)
 * is checked before C's hydrophobic rule, or it could never apply (C is in the hydrophobic group).
 */

import { Color } from "./color.ts";
import { AminoGroupColors as G } from "./dictionaries.ts";

const HYDROPHOBIC = "WLVIMAFCYHP";

/**
 * Jalview's table, a rule per line in order: residues, color, criteria. `KR>60`: K and R together above 60% of the
 * column; `KRQ>85*`: any one of K, R, Q above 85%; `%` is the hydrophobic group.
 */
const RuleTable = `
  AILMFWV hydrophobic %>60
  C       cysteine    C>85
  C       hydrophobic %>60
  KR      positive    KR>60 KRQ>85*
  E       negative    KR>60 QE>50 ED>50 EQD>85*
  D       negative    KR>60 DEN>85* ED>50
  N       polar       N>50 ND>85*
  Q       polar       KR>60 QE>50 QTKR>85*
  ST      polar       %>60 TS>50 ST>85*
  G       glycine     G>0
  P       proline     P>0
  HY      aromatic    %>60 WYACPQFHILMV>85*
`;

interface Criterion {
  readonly letters: string;
  readonly threshold: number;
  /** Each letter's share on its own, not the group's. */
  readonly anyOne: boolean;
}
interface Rule {
  readonly color: Color;
  readonly criteria: readonly Criterion[];
}

/** Rules per residue, in order: the first that holds colors it. */
const Rules: Readonly<Record<string, readonly Rule[]>> = (() => {
  const out: Record<string, Rule[]> = {};
  for (const line of RuleTable.trim().split("\n")) {
    const [residues, group, ...criteria] = line.trim().split(/\s+/);
    const rule: Rule = {
      color: G[group as keyof typeof G],
      criteria: criteria.map((c) => {
        const [letters, threshold] = c.replace("%", HYDROPHOBIC).split(">");
        return { letters, threshold: Number.parseFloat(threshold), anyOne: threshold.endsWith("*") };
      }),
    };
    for (const r of residues) out[r] = [...(out[r] ?? []), rule];
  }
  return out;
})();

/** Letter counts of one column, reused: `counts[code - 65]` for `A`..`Z`, `total` their sum. */
export interface ColumnComposition {
  readonly counts: Int32Array;
  total: number;
}

export function createColumnComposition(): ColumnComposition {
  return { counts: new Int32Array(26), total: 0 };
}

/** Counts the letters of a column (one character per row; gaps and other characters left out) into `out`. */
export function countColumn(column: string, out: ColumnComposition): ColumnComposition {
  out.counts.fill(0);
  out.total = 0;
  for (let i = 0; i < column.length; i++) {
    const code = column.charCodeAt(i) & ~32; // upper case
    if (code >= 65 && code <= 90) {
      out.counts[code - 65]++;
      out.total++;
    }
  }
  return out;
}

function share(c: ColumnComposition, letters: string): number {
  let n = 0;
  for (let i = 0; i < letters.length; i++) n += c.counts[letters.charCodeAt(i) - 65];
  return (100 * n) / c.total;
}

function holds(c: ColumnComposition, { letters, threshold, anyOne }: Criterion): boolean {
  if (!anyOne) return share(c, letters) > threshold;
  for (const letter of letters) if (share(c, letter) > threshold) return true;
  return false;
}

/** The Clustal X color of `letter` in a column of composition `c`; `undefined`: not colored. */
export function clustalXColor(c: ColumnComposition, letter: string): Color | undefined {
  if (c.total === 0) return undefined;
  const rules = Rules[letter.toUpperCase()];
  if (rules === undefined) return undefined;
  for (const { color, criteria } of rules) if (criteria.some((k) => holds(c, k))) return color;
  return undefined;
}

/**
 * The table as `rules` (a named table, `expand`): per criterion a rule on `share` with `in` the line's
 * residues, in the table's order, so a residue's first rule that holds is its first criterion that holds. A group's
 * criterion is one share; "any one of" is a rule per letter. `min` is inclusive, the table's thresholds strict: a
 * share is a count over the column's letters, so `t + 1e-9` is strict for columns of fewer than 10^8 letters.
 */
export function clustalXRules(): unknown[] {
  const out: unknown[] = [];
  for (const line of RuleTable.trim().split("\n")) {
    const [residues, group, ...criteria] = line.trim().split(/\s+/);
    const color = Color.toHex(G[group as keyof typeof G]);
    for (const c of criteria) {
      const [letters, t] = c.replace("%", HYDROPHOBIC).split(">");
      for (const share of t.endsWith("*") ? [...letters] : [letters]) {
        out.push({ when: [{ measure: { share }, min: Number.parseFloat(t) / 100 + 1e-9, in: [...residues] }], color });
      }
    }
  }
  return out;
}
