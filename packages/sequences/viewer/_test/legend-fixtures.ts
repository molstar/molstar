/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Documents for the compact legend (seqspec:interaction#the-legend): many tracks of a few kinds, as track templates
// give them (ST1b-b).

import { createSeqViewSpec } from "@molstar/seq-core/builder";

const A = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ";

/**
 * 25 tracks of three kinds on one sequence, in a template's order: five genes (`gene1`…), ten transcripts by part
 * (`tx1.1`…, UTR and exon) and ten domains (`dom1.1`…). Each track's layer is `l.<track>`. `plain`: the genes have no
 * categories, so nothing to list.
 */
export function kindsDoc(plain = false) {
  const b = createSeqViewSpec({ title: "Kinds" });
  const a = b.data.sequence({ ref: "a", label: "A", value: A });
  const v = b.view({ ref: "v", label: "V" });
  v.block({ sequence: a });
  const track = (ref: string, label: string, start: number[], category?: string[]) => {
    const value = { start, stop: start.map((x) => x + 2), ...(category ? { category } : {}) };
    const data = b.data.annotation({ ref: `d.${ref}`, geometry: "ranges", target: a, value });
    v.track({ ref, label }).ranges({ ref: `l.${ref}`, data });
  };
  for (let g = 1; g <= 5; g++) {
    track(`gene${g}`, `Gene G${g}`, [g], plain ? undefined : ["gene"]);
    for (let k = 1; k <= 2; k++) {
      track(`tx${g}.${k}`, `Transcript G${g}.${k}`, [g, g + 3, g + 6], ["UTR", "exon", "exon"]);
      track(`dom${g}.${k}`, `G${g}.${k} domain`, [g + k], ["domain"]);
    }
  }
  return b.getDocument();
}

/** `n` tracks, each its own category: `n` groups that do not merge (`t0`…, labeled `T0`…). */
export function distinctDoc(n: number) {
  const b = createSeqViewSpec({ title: "Distinct" });
  const a = b.data.sequence({ ref: "a", label: "A", value: A });
  const v = b.view({ ref: "v", label: "V" });
  v.block({ sequence: a });
  for (let k = 0; k < n; k++) {
    const value = { start: [k], stop: [k + 1], category: [`c${k}`] };
    const data = b.data.annotation({ ref: `d${k}`, geometry: "ranges", target: a, value });
    v.track({ ref: `t${k}`, label: `T${k}` }).ranges({ data });
  }
  return b.getDocument();
}
