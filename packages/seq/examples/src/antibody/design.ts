/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Derived } from './derived.js';
import { type Chain, ImgtCdrs, imgtCdrLabels, kabatFromImgt, type Labels } from './numbering.js';

// The design: generated from the wild type by an edit list, never recorded (no public design
// with indels against trastuzumab has a permissive licence). Substitutions in CDR-H2 and CDR-H3, a two-residue
// insertion in CDR-H3 (which gives IMGT's 111.1 and 112.1, out of index order: EC-7) and a one-residue deletion in
// CDR-L1; an edit is placed by the wild type's Kabat label. Its numberings follow the schemes' rules (`numbering.ts`):
// a CDR whose length changed is numbered again by IMGT's rule, the frameworks keep their labels, Kabat's is derived
// from IMGT's as the wild type's is. Its alignment to the wild type is the edit list's, not the numberings': a column
// pairs a residue the design kept with the wild type's, so around an indel the paired residues' labels differ.

/** An edit of the wild type, at a residue by its Kabat label: a substitution, a deletion, or letters inserted after it. */
export type Edit =
  | { readonly chain: Chain; readonly at: string; readonly to: string }
  | { readonly chain: Chain; readonly at: string; readonly delete: true }
  | { readonly chain: Chain; readonly after: string; readonly insert: string };

/** The design's edits (`antibody.design`, an option of the generator): two in CDR-H2, three in CDR-H3, one in CDR-L1. */
export const DesignEdits: readonly Edit[] = [
  // N54Q removes CDR-H2's NG deamidation motif; Y56H
  { chain: 'heavy', at: '54', to: 'Q' },
  { chain: 'heavy', at: '56', to: 'H' },
  { chain: 'heavy', at: '98', to: 'E' },
  { chain: 'heavy', after: '100', insert: 'GY' },
  { chain: 'light', at: '30', delete: true },
];

export interface DesignChain {
  readonly letters: string;
  readonly imgt: Labels;
  readonly kabat: Labels;
  /** The design and the wild type aligned, column by column (`-` a gap): an insertion's columns, a deletion's. */
  readonly aligned: { readonly design: string; readonly wildType: string };
}

export type Design = Readonly<Record<Chain, DesignChain>>;

/** Which IMGT CDR a label is in (0 outside the CDRs). */
function imgtCdr(label: string | null): 0 | 1 | 2 | 3 {
  if (label === null) return 0;
  const n = Number.parseInt(label, 10);
  for (const c of [1, 2, 3] as const) if (n >= ImgtCdrs[c][0] && n <= ImgtCdrs[c][1]) return c;
  return 0;
}

/** The design of a wild type by an edit list. */
export function design(wt: Derived, edits: readonly Edit[] = DesignEdits): Design {
  const one = (chain: Chain): DesignChain => {
    const { letters, imgt } = wt[chain];
    const kabat = kabatFromImgt(chain, imgt);
    // a slot a residue of the wild type: what the design has there (none, its letter, or letters after it)
    const out = Array.from(letters);
    for (const e of edits.filter((e) => e.chain === chain)) {
      const label = 'after' in e ? e.after : e.at;
      const i = kabat.indexOf(label);
      if (i < 0) throw new Error(`${chain}: no residue at Kabat ${label}`);
      if ('to' in e) out[i] = e.to + out[i].slice(1);
      else if ('delete' in e) out[i] = '';
      else out[i] += e.insert;
      // an indel is in a CDR (a framework keeps its labels)
      if (!('to' in e) && imgtCdr(imgt[i]) === 0)
        throw new Error(`${chain}: an indel at Kabat ${label} is not in a CDR`);
    }
    const design: string[] = [];
    const wildType: string[] = [];
    const labels: (string | null)[] = [];
    const cdrs: { cdr: 1 | 2 | 3; at: number[] }[] = [];
    out.forEach((s, i) => {
      const c = imgtCdr(imgt[i]);
      const cdr =
        c === 0 ? undefined : (cdrs.find((r) => r.cdr === c) ?? (cdrs.push({ cdr: c, at: [] }), cdrs[cdrs.length - 1]));
      if (s === '') {
        design.push('-');
        wildType.push(letters[i]);
        return;
      }
      for (let k = 0; k < s.length; k++) {
        cdr?.at.push(labels.length);
        labels.push(imgt[i]);
        design.push(s[k]);
        wildType.push(k === 0 ? letters[i] : '-');
      }
    });
    // each CDR numbered by IMGT's rule for its length; for the wild type's lengths it is the published numbering
    for (const { cdr, at } of cdrs) imgtCdrLabels(cdr, at.length).forEach((l, k) => (labels[at[k]] = l));
    return {
      letters: out.join(''),
      imgt: labels,
      kabat: kabatFromImgt(chain, labels),
      aligned: { design: design.join(''), wildType: wildType.join('') },
    };
  };
  return { light: one('light'), heavy: one('heavy') };
}
