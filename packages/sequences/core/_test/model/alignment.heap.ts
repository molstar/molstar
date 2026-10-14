/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The alignment build at 20k rows (MD1's Opus review): its residues are detected and checked by one counting pass
// over typed arrays, with no string or boxed value per residue. A per-residue allocation over 6M residues would be
// tens of MB; the build's per-alignment objects are a few kB. The build is one run (`heapGrowth`); the check is a
// repeated step (`heapScaling`).

import { describe, expect, it } from "vitest";
import { Alignments, type Issue, Molecules, type ResidueCheck, StringTable } from "../../src/model/index.ts";
import { controlCaught, describeScaling, heapGrowth, heapScaling, isFlat, LIMIT, leak } from "./heap.ts";
import { alignment } from "./support/generators.ts";

const fixture = alignment({ rows: 20_000, length: 300, seed: 5 });
const input = { key: "big", flavor: "aligned" as const, gapped: [...fixture.rows], labels: [...fixture.labels] };

describe("alignment build", () => {
  it("builds a 20k × 300 alignment, detecting its alphabet, without allocating per residue", () => {
    const strings = new StringTable();
    const issues: Issue[] = [];
    let alphabet: string | null = null;
    const build = () => {
      issues.length = 0;
      alphabet = Alignments.build(input, 0, strings, issues, "/alignments/0")?.alphabet ?? null;
    };
    const growth = heapGrowth(build);
    expect(issues).toEqual([]);
    expect(alphabet).toBe("protein");
    expect(growth).toBeLessThan(LIMIT);
    // The template's control: the fault it guards against, an object a residue, is caught by the same measurement.
    const cells = 20_000 * 300;
    expect(
      controlCaught(() =>
        heapGrowth(() => {
          build();
          for (let k = 0; k < cells; k++) leak(k);
        }),
      ),
    ).toBe(true);
  });

  it("checks 6M residues without allocating", () => {
    const strings = new StringTable();
    const a = Alignments.build(input, 0, strings, [], "");
    const residues = a?.rows.residues ?? new Uint8Array(0);
    const check: ResidueCheck = { at: -1, invalid: -1 };
    let alphabet = "";
    const step = () => {
      alphabet = Molecules.checkResidues(residues, undefined, check);
    };
    const r = heapScaling(4, step);
    expect(residues.length).toBeGreaterThan(5_000_000);
    expect([alphabet, check.at, check.invalid]).toEqual(["protein", -1, -1]);
    expect(isFlat(r), describeScaling(r, 4)).toBe(true);
    // The control: an object a residue in the same step is caught (a step is 6M residues: one leak a step is not
    // what this guard is for).
    const perResidue = () => {
      step();
      for (let k = 0; k < residues.length; k++) leak(k);
    };
    expect(controlCaught(() => heapScaling(4, perResidue))).toBe(true);
  });
});
