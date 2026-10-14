/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The edge geometries every renderer is checked on (seqspec:viewer#extension-points; the conformance kit's `edges`),
// as document rewrites a fixture applies to its own document: where draws and hits have gone wrong, a block's extent
// read as the hull of its two ends rather than through its runs (seqspec:model#pieces-and-layout-maps). Sequences and
// alignments given by value are the ones rewritten (an alignment's columns are reversed or filtered only: an axis
// over them is one block, with no origin); a document over
// sources keeps its blocks as they are, and the kit fails a geometry that did not apply unless the fixture skips it.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";

/**
 * - `reversed`: every section's pieces in reverse (`direction: -1`);
 * - `origin`: every sequence cyclic, each unranged block shown from origin 30 (at most half its length), so ranges
 *   cross the shown seam;
 * - `stacked`: every section stacked, its first block followed by a shorter piece of the same sequence;
 * - `filtered`: every section showing two of every three positions of its first block's sequence (a run per kept pair:
 *   a break every third position);
 * - `copies`: every section's first block shown twice side by side, two copies of one unit.
 */
export const EdgeKinds = ["reversed", "origin", "stacked", "filtered", "copies"] as const;
export type EdgeKind = (typeof EdgeKinds)[number];

type Node = { kind: string; ref?: string; params?: Record<string, unknown>; children?: Node[] };
const walk = (n: Node, f: (n: Node) => void) => {
  f(n);
  for (const c of n.children ?? []) walk(c, f);
};

/** `doc` rewritten into one edge geometry (a copy; `doc` is left as it is). */
export function edgeDocument(doc: SeqViewSpecDoc, edge: EdgeKind): SeqViewSpecDoc {
  const d = structuredClone(doc);
  const root = d.root as Node;
  const lengths = new Map<string, number>();
  let data: Node | undefined;
  walk(root, (n) => {
    if (n.kind === "data") data ??= n;
    if (n.kind === "sequence" && n.ref && typeof n.params?.value === "string")
      lengths.set(n.ref, n.params.value.length);
    const rows = n.params?.value;
    if (
      n.kind === "alignment" &&
      n.ref &&
      Array.isArray(rows) &&
      typeof rows[0] === "string" &&
      n.params?.flavor !== "a3m"
    )
      lengths.set(n.ref, rows[0].length);
  });
  const targetOf = (b: Node) => (b.params?.sequence ?? b.params?.alignment) as string;
  const lengthOf = (b: Node) => lengths.get(targetOf(b)) ?? -1;
  walk(root, (n) => {
    if (edge === "origin" && n.kind === "sequence" && n.ref && lengths.has(n.ref))
      n.params = { ...n.params, cyclic: true };
    if (edge === "origin" && n.kind === "block" && lengthOf(n) > 1 && n.params?.range === undefined)
      n.params = { ...n.params, origin: Math.min(30, lengthOf(n) >> 1) };
    if (n.kind !== "section") return;
    if (edge === "reversed") n.params = { ...n.params, direction: -1 };
    const kids = n.children ?? [];
    const at = kids.findIndex((c) => c.kind === "block" && lengthOf(c) > 1);
    if (at < 0) return;
    const first = kids[at];
    const len = lengthOf(first);
    const copy = (params: Record<string, unknown>): Node => ({
      kind: "block",
      ref: `${first.ref ?? n.ref}.edge`,
      params: { ...first.params, ...params },
    });
    // An axis over an alignment is its one block of columns: never stacked, never copied.
    const aligned = first.params?.alignment !== undefined;
    if (edge === "stacked" && !aligned) {
      n.params = { ...n.params, arrangement: "stacked" };
      kids.splice(at + 1, 0, copy({ range: [0, Math.max(1, Math.floor(len / 2))] }));
    }
    if (edge === "copies" && !aligned) {
      // Two copies of one unit are merged into one block unless the section shows each (`copies: "each"`).
      n.params = { ...n.params, copies: "each" };
      kids.splice(at + 1, 0, copy({}));
    }
    if (edge === "filtered" && data) {
      const ref = `${n.ref}.edge.keep`;
      const start = Array.from({ length: Math.ceil(len / 3) }, (_, i) => 3 * i);
      const stop = start.map((s) => Math.min(len, s + 2));
      data.children = [
        ...(data.children ?? []),
        {
          kind: "annotation",
          ref,
          params: { geometry: "ranges", target: targetOf(first), value: { start, stop } },
        },
      ];
      n.params = { ...n.params, positions: ref };
    }
  });
  return d;
}
