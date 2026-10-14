/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The load path at 3M residues (MD1; M2's optional check): a long DNA string builds without a short-lived string or
// array per residue. The old build upper-cased the string (3M characters) and looked each residue up through
// `String.fromCharCode`; the 128-entry tables read the codes straight into the molecule's array. The build is one run
// (`heapGrowth`); detection is a repeated step (`heapScaling`).

import { describe, expect, it } from "vitest";
import { type Issue, Molecules, StringTable } from "../../src/model/index.js";
import { controlCaught, describeScaling, heapGrowth, heapScaling, isFlat, LIMIT, leak } from "./heap.js";
import { prng } from "./support/generators.js";

const SIZE = 3_000_000;
const rand = prng(11);
// Mixed case, so that upper-casing has work to do.
const dna = Array.from({ length: SIZE }, () => "ACGTacgt"[(rand() * 8) | 0]).join("");

describe("molecule build", () => {
  it("builds a 3M-residue DNA string without allocating beyond its own array", () => {
    const strings = new StringTable();
    const issues: Issue[] = [];
    let length = 0;
    const build = () => {
      length = Molecules.build({ key: "chr", residues: dna }, 0, strings, issues, "/molecules/0")?.length ?? 0;
    };
    const growth = heapGrowth(build);
    expect(issues).toEqual([]);
    expect(length).toBe(SIZE);
    expect(growth).toBeLessThan(LIMIT);
    // The template's control: the old build's fault, a short-lived object a residue, is caught by the same measurement.
    const perResidue = () => {
      build();
      for (let k = 0; k < SIZE; k++) leak(k);
    };
    expect(controlCaught(() => heapGrowth(perResidue))).toBe(true);
  });

  it("detects the alphabet of 3M residues without allocating", () => {
    const code = Uint8Array.from(dna.toUpperCase(), (c) => c.charCodeAt(0));
    const found = { at: -1 };
    let alphabet = "";
    const step = () => {
      alphabet = Molecules.detectAlphabet(code, found);
    };
    const r = heapScaling(4, step);
    expect(alphabet).toBe("dna");
    expect(isFlat(r), describeScaling(r, 4)).toBe(true);
    // The control: an object a residue in the same step is caught.
    const perResidue = () => {
      step();
      for (let k = 0; k < SIZE; k++) leak(k);
    };
    expect(controlCaught(() => heapScaling(4, perResidue))).toBe(true);
  });
});
