/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `Systems` (seqspec:model#coordinate-systems; VS-27, GS-13): a system's key and its name for a person are two things
// (an alignment's columns are keyed `msa`, named `msa columns`), a complex is named by its key, a row ref's index is its
// digits only, and a system's numberings are found by name (a namespace alone too) or its default, listed once by name.

import { describe, expect, it } from "vitest";
import { CoordinateSystemKind as K, type KeyTarget, SeqSets, Systems } from "../../src/model/index.js";

function set() {
  return SeqSets.build({
    molecules: [
      {
        key: "p",
        residues: "MKVLA",
        numberings: [{ name: "auth", method: "offset", first: 101, isDefault: true }],
      },
    ],
    units: [
      {
        key: "u",
        molecule: "p",
        numberings: [
          { name: "auth", method: "offset", first: 201 },
          { name: "uniprot:P04637", method: "offset", first: 1 },
        ],
      },
    ],
    complexes: [{ key: "c", units: ["u"] }],
    alignments: [{ key: "msa", gapped: ["MKV-A", "MKVLA"], alphabet: "protein" }],
  } as never).set;
}

describe("Systems", () => {
  it("keys a system as `byName` reads it back, and names it for a person apart", () => {
    const s = set();
    const t: KeyTarget = { kind: -1, id: -1, row: -1 };
    expect(Systems.refOf(s, K.columns, 0, -1)).toBe("msa");
    expect(Systems.nameOf(s, K.columns, 0, -1)).toBe("msa columns");
    expect(Systems.refOf(s, K.row, 0, 1)).toBe("msa/1");
    expect(s.byName(Systems.refOf(s, K.row, 0, 1) as string, t) && t.row).toBe(1);
    expect([Systems.refOf(s, K.complex, 0, -1), Systems.nameOf(s, K.complex, 0, -1)]).toEqual(["c", "c"]);
  });

  it("reads a row ref's index by its digits only", () => {
    expect(["msa/3", "msa/12", "msa/1e1", "msa/0x3", "msa/3.0", "msa/", "msa"].map(Systems.rowOf)).toEqual([
      3, 12, -1, -1, -1, -1, -1,
    ]);
  });

  it("lists a unit's numberings once by name, and finds them by name, by namespace, or the default", () => {
    const s = set();
    expect(Systems.numberingNames(s, K.unit, 0, -1)).toEqual(["auth", "uniprot:P04637", "ordinal"]);
    const first = (names: string[]) => Systems.numberingsOf(s, K.unit, 0, -1, names).map((n) => s.strings.get(n.name));
    expect(first(["uniprot", "auth", "none"])).toEqual(["uniprot:P04637", "auth"]);
    expect(first([])).toEqual(["auth"]);
  });
});
