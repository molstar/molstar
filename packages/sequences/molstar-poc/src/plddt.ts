/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `molstar:plddt` (seqspec:libraries#the-mol-adapter): a predicted model's per-residue confidence as a `values`
// annotation on its complex's polymer units, read from the B-factor column, where AlphaFold DB and ColabFold write
// pLDDT (each residue's atoms carry one value; the first atom's is read). Its default color is AlphaFold DB's four
// bands.

import { type ProviderEntry, RequiredField, ref, SimpleParamsSchema, type ThemeRef } from "@molstar/seq-core/format";
import { chainFor, complexUnits } from "./seq-units.ts";
import { StructureSource, structureOf } from "./source.ts";

/** AlphaFold DB's pLDDT bands (very low < 50 ≤ low < 70 ≤ confident < 90 ≤ very high), as MVS's discrete palette. */
export const PlddtBands: ThemeRef = {
  theme: "field_palette",
  params: {
    field: "value",
    palette: {
      kind: "discrete",
      mode: "absolute",
      colors: [
        ["#FF7D45", 0],
        ["#FFDB13", 50],
        ["#65CBF3", 70],
        ["#0053D6", 90],
      ],
    },
  },
};

export const PlddtProvider: ProviderEntry = {
  name: "molstar:plddt",
  options: SimpleParamsSchema({ of: RequiredField(ref("complex"), "the complex of a predicted model") }),
  geometry: "values",
  defaultColor: PlddtBands,
  async compute(ctx, options) {
    const of = (options as { of: string }).of;
    const parsed = ctx.source(of);
    const src = parsed?.kind === "object" ? parsed.value : null;
    if (!StructureSource.is(src) || !ctx.set) return null;
    const structure = await structureOf(src, 0, null);
    const target: string[] = [];
    const index: number[] = [];
    const value: number[] = [];
    for (const u of complexUnits(ctx.set, of)) {
      const c = chainFor(structure, u.mvs);
      if (c?.kind !== "polymer") continue;
      const b = c.unit.model.atomicConformation.B_iso_or_equiv;
      for (let x = 0; x < c.count; x++) {
        target.push(u.ref);
        index.push(c.labels[x] - u.first);
        value.push(b.value(c.unit.elements[c.start[x]]));
      }
    }
    return { key: "", geometry: "values", semanticType: "plddt", label: "pLDDT", columns: { target, index, value } };
  },
};
