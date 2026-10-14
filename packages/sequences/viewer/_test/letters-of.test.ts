/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Another polymer's letters (`letters` with `of`) follow the letters' rule at every zoom, one size and one form for the
// whole row: through a gapped mapping every mapped residue's letter in view or none; through a codon map every
// residue's three-letter code or none; never a column's text left out among its neighbors' (the dropped columns a fit
// measured per glyph on rounded cell widths gave).

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { recording } from "./composer-fixtures.ts";
import { lit } from "./overlay-fixtures.ts";

const UNIT = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRV";
/** 400 residues; the design lacks 15–17 and has GSGS after 30. */
const WT = UNIT.repeat(10);
const CODING = "ATGGCTTCCAAGGAACTGACCCCACGTGACGGCCATTGGTACTTCAACAGCCAGGTGGAG";
/** 480 coding nucleotides and a stop: 160 residues. */
const GENE = `${CODING.repeat(8)}TAA`;

function mappedDoc(form?: "three") {
  const b = createSeqViewSpec();
  const wt = b.data.sequence({ ref: "wt", alphabet: "protein", value: WT });
  const d = b.data.sequence({
    ref: "d",
    alphabet: "protein",
    value: `${WT.slice(0, 14)}${WT.slice(17, 30)}GSGS${WT.slice(30)}`,
  });
  const gapped: [string, string] = [
    `${WT.slice(0, 30)}----${WT.slice(30)}`,
    `${WT.slice(0, 14)}---${WT.slice(17, 30)}GSGS${WT.slice(30)}`,
  ];
  b.data.mapping({ ref: "d.on_wt", from: wt, to: d, method: "pairwise", gapped });
  b.view({ ref: "v", label: "V" })
    .section({ ref: "v.main" })
    .block({ sequence: wt })
    .track({ label: "" })
    .letters({ of: d, ...(form ? { residue_form: form } : {}) });
  return b.getDocument();
}

function codonDoc() {
  const b = createSeqViewSpec();
  const gene = b.data.sequence({ ref: "gene", alphabet: "dna", value: GENE });
  const protein = translate(b.data, { ref: "protein", from: gene, cds: [[0, GENE.length]] });
  b.view({ ref: "v", label: "V" })
    .section({ ref: "v.main" })
    .block({ sequence: gene })
    .track({ label: "" })
    .letters({ of: protein });
  return b.getDocument();
}

const RESIDUE = /^[A-Z][a-z]{0,2}$/;

/** At each zoom (CSS px a position, from the left end): the residue texts drawn, their centers and fonts. */
async function texts(doc: ReturnType<typeof mappedDoc>, cells: readonly number[]) {
  const t = await lit(doc, { width: 1200, height: 300 });
  const out: { cell: number; strings: string[]; xs: number[]; fonts: Set<number> }[] = [];
  for (const cell of cells) {
    t.vp?.setView(cell, 0, "api");
    const be = t.frame(recording());
    const b = be.textItems;
    const strings: string[] = [];
    const xs: number[] = [];
    for (let i = 0; i < b.n; i++) {
      if (!RESIDUE.test(b.text[i])) continue;
      strings.push(b.text[i]);
      xs.push(b.x[i]);
    }
    const fonts = new Set<number>();
    for (let k = 0; k < be.opCount; k++) {
      const i0 = be.ops[4 * k + 1];
      // A text op's first argument is its font.
      if (be.ops[4 * k + 2] > 0 && i0 < b.n && RESIDUE.test(b.text[i0])) fonts.add(be.args[be.ops[4 * k + 3]]);
    }
    out.push({ cell, strings, xs, fonts });
  }
  return out;
}

/** The steps between consecutive centers, in residues of `span` positions, rounded. */
const steps = (xs: number[], cell: number, span: number) =>
  xs.slice(1).map((x, i) => Math.round((x - xs[i]) / (cell * span)));

const CELLS = [1, 2, 3, 4, 5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 12, 13, 14, 16, 18, 20, 24];

describe("letters of another polymer", () => {
  it("through a gapped mapping: one letter a column, none left out, at one size; or none", async () => {
    let drawn = 0;
    for (const r of await texts(mappedDoc(), CELLS)) {
      if (r.strings.length === 0) continue;
      drawn++;
      expect(
        r.strings.every((s) => s.length === 1),
        `cell ${r.cell}`,
      ).toBe(true);
      expect(r.fonts.size, `cell ${r.cell}`).toBe(1);
      // Every next column, but across the deletion's three.
      expect(
        steps(r.xs, r.cell, 1).filter((x) => x !== 1),
        `cell ${r.cell}`,
      ).toEqual([4]);
    }
    expect(drawn).toBeGreaterThan(8);
  });

  it("through a codon map: every residue's three-letter code, else every one letter, none left out; or none", async () => {
    const forms = new Set<number>();
    for (const r of await texts(codonDoc(), CELLS)) {
      if (r.strings.length === 0) continue;
      const lengths = new Set(r.strings.map((s) => s.length));
      expect(lengths.size, `cell ${r.cell}`).toBe(1);
      forms.add(r.strings[0].length);
      expect(r.fonts.size, `cell ${r.cell}`).toBe(1);
      expect(
        steps(r.xs, r.cell, 3).filter((x) => x !== 1),
        `cell ${r.cell}`,
      ).toEqual([]);
    }
    // Both forms are met across the zooms: one letter where three do not fit a codon's cells.
    expect([...forms].sort()).toEqual([1, 3]);
  });

  it("`residue_form: three` through a gapped mapping: every three-letter code, none left out; or none", async () => {
    let drawn = 0;
    for (const r of await texts(mappedDoc("three"), CELLS)) {
      if (r.strings.length === 0) continue;
      drawn++;
      expect(
        r.strings.every((s) => s.length === 3),
        `cell ${r.cell}`,
      ).toBe(true);
      expect(r.fonts.size).toBe(1);
      expect(
        steps(r.xs, r.cell, 1).filter((x) => x !== 1),
        `cell ${r.cell}`,
      ).toEqual([4]);
    }
    expect(drawn).toBeGreaterThan(2);
  });
});
