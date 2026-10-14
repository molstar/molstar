/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Case 22's triage (the globin family, `core/_test/format/cases.ts`, in the Viewer cases story), shared with the smoke
// test that applies each: the rows layer's `order` and `filter` over its key tables (the review status, `row_stats`,
// `redundancy`), as patches; generic params only.

import { FAMILY_MEMBERS } from "../../../../core/_test/format/cases.js";

export interface TriageMode {
  readonly key: string;
  readonly label: string;
  readonly order: unknown;
  readonly filter: readonly unknown[];
}
const BY_IDENTITY = { data: "stats", field: "identity", direction: "desc" };
/** The layer's `order` and `filter` for each choice; the tree is drawn only in tree order. */
export const TRIAGE_MODES: readonly TriageMode[] = [
  { key: "tree", label: "Tree order", order: "tree", filter: [] },
  {
    key: "reviewed",
    label: "Reviewed first, then by identity",
    order: [{ data: "review", field: "status" }, BY_IDENTITY],
    filter: [],
  },
  {
    key: "reviewedOnly",
    label: "Reviewed only",
    order: BY_IDENTITY,
    filter: [{ data: "review", field: "status", values: ["reviewed"] }],
  },
  {
    key: "covering",
    label: "Covering 80% or more",
    order: BY_IDENTITY,
    filter: [{ data: "stats", field: "coverage", min: 0.8 }],
  },
  {
    key: "representatives",
    label: "One per cluster at 50%",
    order: BY_IDENTITY,
    filter: [{ data: "nr50", field: "representative", values: [true] }],
  },
];

/** The patch that sets a triage mode on the members layer. */
export const triagePatch = (m: TriageMode) => ({
  format: "seq_view_spec_patch" as const,
  ops: [{ op: "set" as const, ref: FAMILY_MEMBERS, params: { order: m.order, filter: m.filter } }],
});
