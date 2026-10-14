/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A wild type and a design through a pairwise mapping, a section each, for the carets' test and heap guard: the design
// has GSGS after WT 20 (its 20–23), WT 30–32 deleted (WT 29–31, the junction before its 33) and a tag of 3 after the
// last (its 41–43); a boundary item at WT 10.

import { createSeqViewSpec } from "@molstar/seq-core/builder";

export const WT = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRV";
/** The design: GSGS after WT 20, WT 30–32 deleted, a tag of 3 after the last. */
export const DESIGN = `${WT.slice(0, 20)}GSGS${WT.slice(20, 29)}${WT.slice(32)}HHH`;

export function caretsDoc() {
  const b = createSeqViewSpec();
  const wt = b.data.sequence({ ref: "wt", value: WT });
  const d = b.data.sequence({ ref: "d", value: DESIGN });
  b.data.mapping({
    ref: "d.on_wt",
    from: wt,
    to: d,
    method: "pairwise",
    gapped: [
      `${WT.slice(0, 20)}----${WT.slice(20)}---`,
      `${WT.slice(0, 20)}GSGS${WT.slice(20, 29)}---${WT.slice(32)}HHH`,
    ],
  });
  const marks = b.data.annotation({ ref: "marks", geometry: "ranges", target: wt, value: { boundary: [10] } });
  const v = b.view({ ref: "v", label: "V" });
  v.section({ ref: "v.wt" })
    .block({ sequence: wt })
    .track({ ref: "t-wt", label: "WT" })
    .letters()
    .ranges({ data: marks });
  v.section({ ref: "v.d" }).block({ sequence: d }).track({ ref: "t-d", label: "D" }).letters();
  return b.getDocument();
}
