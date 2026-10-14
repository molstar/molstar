/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
// Tables with insertion codes (seqspec:format#positions): `ins_code` (and `b_ins_code`) completes an author number in
// `pos` or `beg` (`b_beg`) of one residue, "", "?" and "." none; a file's columns renamed onto them by `fields`.
// Without codes a table reads as before. The errors say what to write instead.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { SeqSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { compile } from "../../src/format/index.js";
import { fixtureLoaders, registries } from "./load-fixtures.js";

const SEQ = "WVARIYPTNGYT";
/** Author numbers 50, 51, 52, 52A, 52B, 53, 54, 55, 56, 57, 58, 59. */
const AUTH: (number | string)[] = [50, 51, 52, "52A", "52B", 53, 54, 55, 56, 57, 58, 59];
const at = (l: string) => AUTH.findIndex((x) => String(x) === l);
const k = "k";

async function compiled(
  geometry: string,
  table: Record<string, unknown> | null,
  file?: { text: string; fields: Record<string, string> },
) {
  const b = createSeqViewSpec();
  const h = b.data.sequence({ ref: "h", value: SEQ }).numbering({ name: "auth", method: "labels", labels: AUTH });
  b.data.sequence({ ref: "k", value: SEQ }).numbering({ name: "auth", method: "labels", labels: AUTH });
  const pairs = geometry === "pairs" ? { b_numbering: "auth" } : {};
  b.data.annotation({
    ref: "a",
    geometry,
    target: h,
    numbering: "auth",
    ...(file ? { uri: "https://example.org/a.tsv", format: "tsv", fields: file.fields } : { value: table }),
    ...pairs,
  } as never);
  b.view({ ref: "v", label: "V" }).block({ sequence: h });
  const r = await compile(
    b.getDocument(),
    registries(false),
    fixtureLoaders({ "https://example.org/a.tsv": file?.text ?? "" }),
  );
  const set = r.set as SeqSet;
  const a = set.annotations.find((x) => x && set.strings.get(x.key) === "a");
  return {
    issues: r.issues.map((i) => [i.code, i.severity, i.message]),
    start: Array.from(a?.start ?? []),
    stop: Array.from(a?.stop ?? []),
    bStart: Array.from(a?.pairs?.bStart ?? []),
  };
}

describe("insertion codes in tables", () => {
  it("without codes a table reads as before; by number alone, 52 is 52", async () => {
    const r = await compiled("values", { pos: [52, 53], value: [1, 2] });
    expect(r).toMatchObject({ issues: [], start: [at("52"), at("53")] });
  });

  it("ins_code completes pos; empty, ? and . are none", async () => {
    const r = await compiled("values", {
      pos: [52, 52, 52, 53, 54, 55],
      ins_code: [null, "A", "B", "?", ".", ""],
      value: [1, 2, 3, 4, 5, 6],
    });
    expect(r).toMatchObject({ issues: [], start: ["52", "52A", "52B", "53", "54", "55"].map(at) });
  });

  it("on a range's beg alone, or beg and end the same, the code names that one residue", async () => {
    const r = await compiled("ranges", { beg: [52, 52], end: [null, 52], ins_code: ["A", "B"] });
    expect(r).toMatchObject({ issues: [], start: [at("52A"), at("52B")], stop: [at("52A") + 1, at("52B") + 1] });
  });

  it("b_ins_code completes a pair's other end", async () => {
    const r = await compiled("pairs", {
      beg: [52],
      ins_code: ["A"],
      b_target: [k],
      b_beg: [52],
      b_ins_code: ["B"],
    });
    expect(r).toMatchObject({ issues: [], start: [at("52A")], bStart: [at("52B")] });
  });

  it("a file's columns renamed by fields: auth_seq_id and pdbx_PDB_ins_code", async () => {
    const text = "auth_seq_id\tpdbx_PDB_ins_code\tscore\n52\t?\t1\n52\tA\t2\n52\tB\t3\n53\t?\t4\n";
    const r = await compiled("values", null, {
      text,
      fields: { pos: "auth_seq_id", ins_code: "pdbx_PDB_ins_code", value: "score" },
    });
    expect(r).toMatchObject({ issues: [], start: ["52", "52A", "52B", "53"].map(at) });
  });

  it("a code that names no residue leaves its row out, saying what the labels are", async () => {
    const r = await compiled("values", { pos: [52, 52], ins_code: ["A", "C"], value: [1, 2] });
    expect(r.start).toEqual([at("52A")]);
    expect(r.issues).toEqual([
      [
        "svs.position.unknown",
        "warning",
        'Annotation "a", row 1: no position "52C" in numbering "auth" (its labels run 50-59; on 52: 52, 52A, 52B).',
      ],
    ]);
  });

  it("instructive errors: a code on a range, on a label that has one, not letters", async () => {
    const r = await compiled("ranges", {
      beg: [50, "52A", 52, 53],
      end: [55, null, null, null],
      ins_code: ["A", "B", "1", "A"],
    });
    expect(r.issues).toEqual([
      [
        "svs.position.form",
        "error",
        'Annotation "a", row 0: ins_code "A" names one residue, but beg 50 to end 55 is a range: write a range\'s ends as whole labels (beg "50A"), or give beg alone.',
      ],
      [
        "svs.position.form",
        "error",
        'Annotation "a", row 1: beg "52A" is not an author number to complete with ins_code "B": give the number alone (beg "52", ins_code "A"), or the whole label without ins_code.',
      ],
      [
        "svs.position.syntax",
        "error",
        'Annotation "a", row 2: insertion code "1" is not letters ("A", "AA"); none is "", "?" or ".".',
      ],
      [
        "svs.position.unknown",
        "warning",
        'Annotation "a", row 3: no position "53A" in numbering "auth" (its labels run 50-59; on 53: 53).',
      ],
    ]);
  });

  it("ins_code with indices is a form error that says to write labels", async () => {
    const r = await compiled("values", { index: [2], ins_code: ["A"], value: [1] });
    expect(r.issues[0]?.slice(0, 2)).toEqual(["svs.position.form", "error"]);
    expect(r.issues[0]?.[2]).toMatch(/ins_code completes an author number in beg or pos .* give positions as labels/);
  });
});

describe("a label range written backwards, or naming a label that is not there (CF-01)", () => {
  it("both ends placed, the end first: svs.position.range, saying to write beg ≤ end", async () => {
    for (const [beg, end] of [
      [55, 53],
      ["52A", 52],
      ["52C", "52A"], // 52C is absent, placed by order at 53: still after 52A
    ]) {
      const r = await compiled("ranges", { beg: [beg], end: [end] });
      expect(r.issues, `${beg}-${end}`).toEqual([
        [
          "svs.position.range",
          "error",
          `Annotation "a", row 0: beg "${beg}" comes after end "${end}" in numbering "auth": write beg ≤ end (a range wraps only on a cyclic sequence).`,
        ],
      ]);
    }
  });

  /** The issues of one range by `labels` (`cyclic`: on a cyclic sequence). */
  const ranged = async (labels: (number | string)[], beg: unknown, end: unknown, cyclic = false) => {
    const b = createSeqViewSpec();
    const value = "WVARIYPTNG".slice(0, labels.length);
    const h = b.data
      .sequence({ ref: "h", value, ...(cyclic ? { cyclic: true } : {}) } as never)
      .numbering({ name: "auth", method: "labels", labels });
    b.data.annotation({
      ref: "a",
      geometry: "ranges",
      target: h,
      numbering: "auth",
      value: { beg: [beg], end: [end] },
    });
    b.view({ ref: "v", label: "V" }).block({ sequence: h });
    const r = await compile(b.getDocument(), registries(false), fixtureLoaders({}));
    return r.issues.map((i) => `${i.code} ${i.message}`);
  };

  it("labels not in index order (IMGT's insertions, a circular permutant) are ordered by index (LA-06)", async () => {
    const imgt = [110, 111, "111.1", "112.1", 112, 113];
    expect(await ranged(imgt, "112.1", 112)).toEqual([]);
    expect(await ranged(imgt, 112, "112.1")).toEqual([
      expect.stringMatching(/^svs\.position\.range .*beg "112" comes after end "112\.1" .*write beg ≤ end/),
    ]);
    const permutant = [50, 51, 52, 53, 1, 2, 3];
    expect(await ranged(permutant, 51, 2)).toEqual([]);
    expect(await ranged(permutant, 2, 51)).toEqual([
      expect.stringMatching(/^svs\.position\.range .*beg "2" comes after end "51"/),
    ]);
  });

  it("on a cyclic sequence, a range across the origin is less than a whole lap, and says so (LA-20)", async () => {
    const ring = [1, 2, 3, 4, 5, 6];
    expect(await ranged(ring, 4, 2, true)).toEqual([]);
    expect(await ranged(ring, 4, 3, true)).toEqual([
      expect.stringMatching(/^svs\.position\.range .*beg "4" .* end "3" .*less than a whole lap/),
    ]);
  });

  it("a label that is not there is named, beg or end, with the numbering's labels", async () => {
    const after = await compiled("ranges", { beg: [70], end: [75] });
    expect(after.issues).toEqual([
      ["svs.position.unknown", "warning", expect.stringMatching(/no position "70" in numbering "auth"/)],
    ]);
    const before = await compiled("ranges", { beg: [40], end: [45] });
    expect(before.issues).toEqual([
      ["svs.position.unknown", "warning", expect.stringMatching(/no position "45" in numbering "auth"/)],
    ]);
  });
});
