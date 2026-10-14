/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `grid` geometry and the `dot_plot` layer type (seqspec:model#annotations, seqspec:format#layer-types): a value
// per pair of positions, a-major, its b side placed as a pair's b end; a dot plot over two blocks of one section.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { KeyKind, type KeyTarget } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { compile, resolveView, supportedAll, validate } from "../../src/format/index.ts";
import { case8 } from "./cases.ts";
import { files, fixtureLoaders, registries } from "./load-fixtures.ts";

const reg = registries();
const t: KeyTarget = { kind: -1, id: -1, row: -1 };

/** Two sequences, s (4) and q (3), and a grid over them: s[0, 2) × q whole, and q × s whole (one missing value). */
function doc(value?: Record<string, unknown>) {
  const b = createSeqViewSpec({ title: "grid" });
  const s = b.data.sequence({ ref: "s", value: "MKVL" });
  const q = b.data.sequence({ ref: "q", value: "AGI" });
  const g = b.data.annotation({
    ref: "g",
    geometry: "grid",
    value: value ?? {
      target: ["s", "q"],
      start: [0, null],
      stop: [2, null],
      b_target: ["q", "s"],
      value: [
        [1, 2, 3, 4, 5, 6],
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, null],
      ],
    },
  });
  const v = b.view({ ref: "v", label: "V" });
  const sec = v.section({ ref: "sec" });
  sec.block({ sequence: s });
  sec.block({ sequence: q });
  sec.track({ label: "Plot" }).dot_plot({ data: g });
  return b.getDocument();
}

describe("grid annotations", () => {
  it("hold a value per pair of positions, a-major; a side with no range is its target's whole domain", async () => {
    const { set, issues } = await compile(doc(), reg, fixtureLoaders());
    expect(issues).toEqual([]);
    expect(set.byName("g", t) && t.kind === KeyKind.annotation).toBe(true);
    const a = set.annotations[t.id];
    expect(a?.geometry).toBe("grid");
    expect([Array.from(a?.start ?? []), Array.from(a?.stop ?? [])]).toEqual([
      [0, 0],
      [2, 3],
    ]);
    expect([Array.from(a?.pairs?.bStart ?? []), Array.from(a?.pairs?.bStop ?? [])]).toEqual([
      [0, 0],
      [3, 4],
    ]);
    expect(a?.grid?.rows.map((r) => r.length)).toEqual([6, 12]);
    // Row 0: s[1] × q[2] is the second a row's third value.
    expect(a?.grid?.rows[0][1 * 3 + 2]).toBe(6);
    expect(Number.isNaN(a?.grid?.rows[1][11] ?? 0)).toBe(true);
    expect(a?.sparse).toBe(null);
  });

  it("leave out a row whose values do not fill its block, and say how many are needed", async () => {
    const { issues } = await compile(
      doc({ target: ["s"], b_target: ["q"], start: [0], stop: [2], value: [[1, 2, 3]] }),
      reg,
      fixtureLoaders(),
    );
    expect(issues.map((i) => [i.code, i.message])).toEqual([
      [
        "svs.columns.type",
        "Annotation \"g\", row 0: 2 × 3 = 6 values are needed (side a's positions × side b's), a-major; got 3.",
      ],
    ]);
  });

  it("keep a row's Float32Array as given, and say a grid's sides never wrap a cyclic origin", async () => {
    const cells = Float32Array.from({ length: 6 }, (_, k) => k / 10);
    // In the document as given (the builder writes JSON values).
    const d = doc({ target: ["s"], b_target: ["q"], start: [0], stop: [2], value: [[]] });
    const g = (d.root.children ?? []).flatMap((n) => n.children ?? []).find((n) => n.ref === "g");
    if (!g?.params) throw new Error("no grid");
    (g.params.value as { value: unknown[] }).value = [cells];
    const { set, issues } = await compile(d, reg, fixtureLoaders());
    expect(issues).toEqual([]);
    expect(set.byName("g", t) && set.annotations[t.id]?.grid?.rows[0]).toBe(cells);
    const b = createSeqViewSpec({ title: "cyclic" });
    b.data.sequence({ ref: "c", value: "ACGTACGTAC", alphabet: "dna", cyclic: true });
    b.data.annotation({ ref: "g", geometry: "grid", value: { target: ["c"], start: [8], stop: [2], value: [[]] } });
    b.view({ ref: "v", label: "V" }).track({ label: "P" }).dot_plot({ data: "g" });
    const wrapped = await compile(b.getDocument(), reg, fixtureLoaders());
    expect(wrapped.issues.map((i) => i.message)).toEqual([
      expect.stringContaining(
        "a grid row's sides never wrap a cyclic origin: write the run across the origin as two rows, [8, 10) and [0, 2)",
      ),
    ]);
  });

  it("an aspect ratio is above 0, and says what to write", () => {
    const d = doc();
    const layer = JSON.stringify(d).replace('"type":"dot_plot"', '"type":"dot_plot","aspect_ratio":0');
    expect(validate(JSON.parse(layer), reg).map((i) => [i.code, i.message])).toEqual([
      [
        "svs.params.type",
        "aspect_ratio 0 is not above 0: write the plot's width ÷ height (1: square, 2: twice as wide), or null for the track's own height",
      ],
    ]);
  });

  it("check a value cell as a list of numbers statically", () => {
    const issues = validate(doc({ target: ["s"], b_target: ["q"], value: [3] }), reg);
    expect(issues.map((i) => i.message)).toEqual([
      expect.stringContaining("a list of float | null (or a Float32Array, a Float64Array), one per pair of positions"),
    ]);
  });
});

describe("dot plots", () => {
  it("place both ends on the section's blocks and color a grid by the numeric default over its domain", async () => {
    const d = doc();
    const { set } = await compile(d, reg, fixtureLoaders());
    const { view, issues } = resolveView(d, set, "v", reg, supportedAll(reg));
    expect(issues).toEqual([]);
    // Side a's two slots and side b's two ends, each on its own block.
    const pl = view.placements;
    const placed = Array.from({ length: pl.count }, (_, i) => [pl.block[i], pl.slot[i]]).sort();
    expect(placed).toEqual([
      [0, 0],
      [0, 3],
      [1, 1],
      [1, 2],
    ]);
    expect(view.layers.themes[0][0]).toMatchObject({
      theme: "field_palette",
      params: { palette: { kind: "continuous", colors: "Blues", value_domain: [0, 10] } },
    });
  });

  it("take pairs or grids, nothing else", () => {
    const b = createSeqViewSpec({ title: "x" });
    const s = b.data.sequence({ ref: "s", value: "MKVL" });
    const r = b.data.annotation({ ref: "r", geometry: "ranges", target: s, value: { start: [0], stop: [2] } });
    b.view({ ref: "v", label: "V" }).track({ label: "Plot" }).dot_plot({ data: r });
    expect(validate(b.getDocument(), reg).map((i) => [i.code, i.message])).toEqual([
      ["svs.layer.data", "dot_plot takes pairs, grid; got ranges"],
    ]);
  });
});

describe("case 8's predicted aligned error", () => {
  it("reads AlphaFold DB's PAE into a grid per pair of chains, scored across, aligned on down; the plot on both", async () => {
    const d = case8();
    const { set, issues } = await compile(d, reg, fixtureLoaders());
    expect(issues.filter((i) => i.severity === "error")).toEqual([]);
    expect(set.byName("pae", t) && t.kind === KeyKind.annotation).toBe(true);
    const a = set.annotations[t.id];
    const key = (u: number) => set.strings.get(set.units.key[u]);
    expect(Array.from(a?.slot ?? [], (s) => key(a?.slots.unit[s] ?? -1))).toEqual([
      "model/A",
      "model/A",
      "model/B",
      "model/B",
    ]);
    expect(Array.from(a?.pairs?.bTarget.unit ?? [], key)).toEqual(["model/A", "model/B", "model/A", "model/B"]);
    // The fixture's row 0 (aligned on barnase 1) at column 1 (scored barnase 2): side a 1, side b 0 of A × A.
    const pae = JSON.parse(
      files()["https://alphafold.ebi.ac.uk/files/AF-BARNASE-BARSTAR-predicted_aligned_error_v4.json"],
    );
    const rows = pae[0].predicted_aligned_error as number[][];
    const nb = (a?.pairs?.bStop[0] ?? 0) - (a?.pairs?.bStart[0] ?? 0);
    expect(a?.grid?.rows[0][1 * nb + 0]).toBe(Math.fround(rows[0][1]));
    const { view } = resolveView(d, set, "chains", reg, supportedAll(reg));
    const k = Array.from(view.layers.type).findIndex((x) => x >= 0 && set.strings.get(x) === "dot_plot");
    expect(k).toBeGreaterThanOrEqual(0);
    const pl = view.placements;
    const blocks = new Set<number>();
    for (let i = 0; i < pl.count; i++) if (pl.layer[i] === k) blocks.add(pl.block[i]);
    expect(blocks.size).toBe(2);
  });
});
