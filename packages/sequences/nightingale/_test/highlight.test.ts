/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// CP3's highlight on codons-fixtures' gene (case 14's, with a protein section and a chromosome section): a manager's
// `highlight` is the axis conversion of the shared projection's runs, from the same loci a host passes (`element`, and
// `containing` for a pointer hover). The old stack's was compared until SW1 (its parity layer, deleted).

import { OrderedSet, type UnitId } from "@molstar/seq-core/model";
import { beforeAll, describe, expect, it } from "vitest";
import { gene, open } from "../../core/_test/format/project-fixtures.ts";
import { axesOf } from "../src/axis.ts";
import type { AxisMeta } from "../src/index.ts";
import { highlightOf } from "./highlight-of.ts";

let g: Awaited<ReturnType<typeof open>>;
let managers: AxisMeta[];
beforeAll(async () => {
  g = await open(gene({ chromosome: true }));
  managers = Array.from({ length: g.view.sections.count }, (_, s) => axesOf(g.view, g.set, s)).flat();
});
const seq = (key: string, a: number, b: number) => ({
  kind: "seq-loci" as const,
  set: g.set,
  elements: [{ unit: g.system(key).id as UnitId, indices: OrderedSet.ofRange(a, b) as never }],
});
/** Every manager's highlight for loci. */
const lit = (loci: ReturnType<typeof seq> | null, closure?: "containing") =>
  managers.map((m) => highlightOf(g.p, m, loci, closure));

describe("a host's loci on each manager", () => {
  it("a residue lights its codon on the gene block, split by the intron, and across the link on the chromosome", () => {
    // Residue 11 (H12) is nucleotides 45 and 76..77 (0-based) of the gene, 100 further on the chromosome.
    expect(lit(seq("protein", 11, 12))).toEqual(["46:46,77:78", "12:12", "146:146,177:178"]);
    expect(lit(seq("protein", 0, 1))).toEqual(["13:15", "1:1", "113:115"]);
  });

  it("a nucleotide lights itself and its residue; with the containing closure, its whole codon too", () => {
    expect(lit(seq("gene/unit", 77, 78))).toEqual(["78:78", "12:12", "178:178"]);
    // A pointer hover (`containing`): the residue that holds it is closed in, so every manager shows the codon.
    const codon = ["46:46,77:78", "12:12", "146:146,177:178"];
    expect(lit(seq("gene/unit", 77, 78), "containing")).toEqual(codon);
    expect(lit(seq("chr/unit", 177, 178), "containing")).toEqual(codon);
  });

  it("nothing is lit by no loci, or by loci of another set", () => {
    expect(lit(null)).toEqual(["", "", ""]);
    const other = { ...seq("protein", 0, 1), set: { ...g.set } };
    expect(lit(other as never)).toEqual(["", "", ""]);
  });
});
