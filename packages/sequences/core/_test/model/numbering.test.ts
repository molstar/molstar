/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from "vitest";
import {
  type Issue,
  type NameId,
  type NameNumbering,
  NO_LABEL,
  Numbering,
  type NumberingInput,
  type PositionLabel,
  type SeqIndex,
  StringTable,
} from "../../src/model/index.js";

const strings = new StringTable();

function numbering(input: NumberingInput, length: number) {
  const issues: Issue[] = [];
  const n = Numbering.build(input, length, strings, issues, "/n");
  expect(issues).toEqual([]);
  return n as Numbering;
}
const labels = (xs: (number | string | null)[]) => numbering({ name: "l", method: "labels", labels: xs }, xs.length);

/** A position label from text, as the compiler writes it: "27A" → (27, "A"). */
function pos(label: number | string): PositionLabel {
  const m = /^(-?\d+)(.*)$/.exec(String(label)) as RegExpExecArray;
  return { num: Number(m[1]), suffix: m[2] ? strings.intern(m[2]) : (-1 as NameId), name: -1 as NameId };
}
function indexOf(n: Numbering, label: number | string) {
  const p = pos(label);
  return Numbering.indexOf(n, p.num, p.suffix);
}
function rangeOf(n: Numbering, beg: number | string, end: number | string): [number, number] | null {
  const out = new Int32Array(2);
  return Numbering.rangeOf(n, pos(beg), pos(end), out) ? [out[0], out[1]] : null;
}
const formatAll = (n: Numbering) =>
  Array.from({ length: n.length }, (_, i) => Numbering.format(n, i as SeqIndex, strings));

describe("numberings", () => {
  // The two levels (a unit's numbering wins over its molecule's: `SeqSets.numberingOf`) are the set's: set.test with MD2.

  it("finds insertion codes and formats them", () => {
    const n = labels([1, 2, "2A", "2B", 3, null, 5]);
    expect(formatAll(n)).toEqual(["1", "2", "2A", "2B", "3", "", "5"]);
    expect([indexOf(n, 2), indexOf(n, "2A"), indexOf(n, "2B"), indexOf(n, "2C"), indexOf(n, 4)]).toEqual([
      1, 2, 3, -1, -1,
    ]);
    const out: PositionLabel = { num: 0, suffix: -1 as NameId, name: -1 as NameId };
    expect(Numbering.labelOf(n, 2 as SeqIndex, out)).toBe(true);
    expect([out.num, strings.get(out.suffix)]).toEqual([2, "A"]);
    expect(Numbering.labelOf(n, 5 as SeqIndex, out)).toBe(false);
    expect(out.num).toBe(NO_LABEL);
  });

  it("orders IMGT's 111.1 and 112B < 112A < 112 by index, and places absent endpoints by comparison", () => {
    const n = labels([110, 111, "111.1", "111.2", "112B", "112A", 112, 113]);
    expect(indexOf(n, "111.1")).toBe(2);
    expect(rangeOf(n, "111.1", "112A")).toEqual([2, 6]);
    expect(rangeOf(n, "112B", 112)).toEqual([4, 7]);
    expect(rangeOf(n, 112, "112B")).toBeNull();
    expect(rangeOf(n, "111.3", 113)).toEqual([4, 8]);
    expect(rangeOf(n, 110, "111.5")).toEqual([0, 4]);
  });

  it("takes insertion codes of several letters (past Z), ordered A … Z, AA …; absent ones placed by that order", () => {
    const n = labels([100, "100A", "100Z", "100AA", "100AB", 101]);
    expect(formatAll(n)).toEqual(["100", "100A", "100Z", "100AA", "100AB", "101"]);
    expect([indexOf(n, "100AA"), indexOf(n, "100AB")]).toEqual([3, 4]);
    expect(rangeOf(n, "100B", "100AA")).toEqual([2, 4]);
    expect(rangeOf(n, "100Z", "100AC")).toEqual([2, 5]);
  });

  it("ends Kabat's CDR-H1 31–35B at 35 when 35A and 35B are absent", () => {
    expect(rangeOf(labels([29, 30, 31, 32, 33, 34, 35, 36]), 31, "35B")).toEqual([2, 7]);
    expect(rangeOf(labels([30, 31, 32, 33, 34, 35, "35A", "35B", 36]), 31, "35B")).toEqual([1, 8]);
    expect(rangeOf(labels([29, 30, 36, 37]), 31, "35B")).toBeNull();
  });

  it("numbers a precursor's mature chain from 1, skipping zero", () => {
    const n = numbering({ name: "mature", method: "offset", first: -47, skipZero: true }, 157);
    expect([0, 46, 47, 156].map((i) => Numbering.format(n, i as SeqIndex, strings))).toEqual(["-47", "-1", "1", "110"]);
    expect([indexOf(n, -47), indexOf(n, -1), indexOf(n, 0), indexOf(n, 1), indexOf(n, 110), indexOf(n, 111)]).toEqual([
      0, 46, -1, 47, 156, -1,
    ]);
    expect(indexOf(n, "1A")).toBe(-1);
    expect(rangeOf(n, -2, 2)).toEqual([45, 49]);
    expect(rangeOf(n, 0, 0)).toBeNull();
  });

  it("maps segments with holes both ways", () => {
    const n = numbering(
      {
        name: "uniprot:X",
        method: "segments",
        segments: [
          [0, 100, 10],
          [15, 200, 5],
        ],
      },
      20,
    );
    expect([
      indexOf(n, 100),
      indexOf(n, 109),
      indexOf(n, 110),
      indexOf(n, 204),
      indexOf(n, 205),
      indexOf(n, 99),
    ]).toEqual([0, 9, -1, 19, -1, -1]);
    expect(formatAll(n).slice(8, 17)).toEqual(["108", "109", "", "", "", "", "", "200", "201"]);
    expect(rangeOf(n, 105, 202)).toEqual([5, 18]);
  });

  it("finds a label in a segment that a later one does not cover, whatever the order of the labels", () => {
    // Index order is not label order: the second segment's labels come first.
    const n = numbering(
      {
        name: "s",
        method: "segments",
        segments: [
          [0, 100, 10],
          [10, 1, 10],
          [25, 300, 5],
        ],
      },
      30,
    );
    expect([indexOf(n, 100), indexOf(n, 109), indexOf(n, 1), indexOf(n, 10), indexOf(n, 11), indexOf(n, 99)]).toEqual([
      0, 9, 10, 19, -1, -1,
    ]);
    expect([indexOf(n, 300), indexOf(n, 304), indexOf(n, 305)]).toEqual([25, 29, -1]);
    // Labels that touch are disjoint.
    const touching = numbering(
      {
        name: "t",
        method: "segments",
        segments: [
          [0, 1, 300],
          [300, 301, 50],
        ],
      },
      350,
    );
    expect([indexOf(touching, 200), indexOf(touching, 301)]).toEqual([199, 300]);
  });

  it("rejects segments whose label ranges overlap or nest", () => {
    // Labels 1..300 and 100..110: the lookup of 200 used to answer "absent" from the nested segment.
    for (const segments of [
      [
        [0, 1, 300],
        [300, 100, 11],
      ],
      [
        [0, 1, 10],
        [10, 10, 10],
      ],
      [
        [0, 50, 10],
        [10, 1, 51],
      ],
    ] as [number, number, number][][]) {
      const issues: Issue[] = [];
      expect(Numbering.build({ name: "s", method: "segments", segments }, 400, strings, issues, "/n")).toBeNull();
      expect(issues.map((i) => [i.code, i.path])).toEqual([["model.numbering.invalid", "/n"]]);
      expect(issues[0].message).toContain("label");
    }
  });

  it("refuses a label outside int32, never wrapped or read as unlabeled (CM-14)", () => {
    for (const bad of [3_000_000_000, -2147483648, "3000000000A"]) {
      const issues: Issue[] = [];
      const built = Numbering.build({ name: "l", method: "labels", labels: [1, bad, 3] }, 3, strings, issues, "/n");
      expect([built, issues.map((i) => i.code), issues[0]?.message]).toEqual([
        null,
        ["model.numbering.invalid"],
        expect.stringContaining("outside the labels a numbering holds"),
      ]);
    }
  });

  it("names the first position of a repeated label", () => {
    const n = labels([5, 6, 5, "5A", 6, "5A"]);
    expect([indexOf(n, 5), indexOf(n, 6), indexOf(n, "5A"), indexOf(n, 7)]).toEqual([0, 1, 3, -1]);
  });

  it("looks names up, and ranges between names", () => {
    const n = numbering({ name: "atom", method: "names", names: ["FE", "NA", null, "NB"] }, 4) as NameNumbering;
    expect(Numbering.indexOfName(n, strings.find("NB"))).toBe(3);
    expect(Numbering.indexOfName(n, strings.intern("XX"))).toBe(-1);
    const name = (s: string): PositionLabel => ({ num: NO_LABEL, suffix: -1 as NameId, name: strings.find(s) });
    const out = new Int32Array(2);
    expect(Numbering.rangeOf(n, name("NA"), name("NB"), out) && [...out]).toEqual([1, 4]);
    expect(formatAll(n)).toEqual(["FE", "NA", "", "NB"]);
  });

  it("uses a sorted table when labels are sparse", () => {
    const n = labels([-5000, 1, 100000, "100000A"]);
    expect((n as { lookupKeys: Int32Array | null }).lookupKeys).not.toBeNull();
    expect([indexOf(n, -5000), indexOf(n, 100000), indexOf(n, "100000A"), indexOf(n, 2)]).toEqual([0, 2, 3, -1]);
  });

  it("reports bad numberings and leaves them out", () => {
    const issues: Issue[] = [];
    const set = Numbering.buildSet(
      [
        { name: "short", method: "labels", labels: [1, 2] },
        { name: "bad", method: "labels", labels: [1, "2A1", 3] },
        { name: "ordinal", method: "offset", first: 0 },
        {
          name: "overlap",
          method: "segments",
          segments: [
            [0, 1, 2],
            [1, 5, 1],
          ],
        },
        { name: "ok", method: "labels", labels: [1, 2, 3] },
      ],
      3,
      strings,
      issues,
      "/numberings",
      true,
    );
    expect(issues.map((i) => [i.code, i.path])).toEqual([
      ["model.numbering.invalid", "/numberings/0"],
      ["model.numbering.invalid", "/numberings/1"],
      ["model.numbering.name", "/numberings/2"],
      ["model.numbering.invalid", "/numberings/3"],
    ]);
    expect(set.items.map((n) => strings.get(n.name))).toEqual(["ordinal", "ok"]);
  });
});
