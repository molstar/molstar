/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for resolve's per-row and per-position paths (VW4), on RB0b's stepped `heapScaling` (a step is one row):
// `cover` (what `positions` and an annotation focus mark) through a codon map, the filter's test of a row's value
// (a number column compared as numbers, no double boxed), and a row's color through the layer's one categorical
// palette. Resolve allocates per layer and per call, never per row or position.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { SeqSet } from "@molstar/seq-core/model";
import { beforeAll, describe, expect, it } from "vitest";
import { compile, type End, layerPalette, paletteValue, resolveView, supportedAll } from "../../src/format/index.ts";
import { cover, matcher } from "../../src/format/resolve/reach.ts";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "../model/heap.ts";
import { fixtureLoaders, registries } from "./load-fixtures.ts";

/** Steps in the shorter run; the longer runs four times as many. */
const N = 20_000;
const ROWS = 4 * N;
const reg = registries();
let set: SeqSet;
let view: ReturnType<typeof resolveView>["view"];

beforeAll(async () => {
  // A 300-residue protein over its 900-nucleotide gene; ROWS one-residue rows on it, with categories and values.
  const b = createSeqViewSpec();
  const protein = Array.from({ length: 300 }, (_, i) => "ACDEFGHIKLMNPQRSTVWY"[i % 20]).join("");
  const gene = b.data.sequence({ ref: "gene", alphabet: "dna", value: "ATG".repeat(300) });
  b.data.sequence({ ref: "p", alphabet: "protein", value: protein });
  b.data.mapping({ ref: "m", from: "p", to: "gene", method: "codon", cds: [[0, 900]] } as never);
  const start = Array.from({ length: ROWS }, (_, i) => i % 300);
  const category = start.map((i) => ["helix", "strand", "turn"][i % 3]);
  const value = start.map((i) => i / 300);
  b.data.annotation({ ref: "a", geometry: "ranges", target: "p", value: { start, category, value } });
  b.view({ ref: "v", label: "V" }).block({ sequence: gene }).track({ label: "A" }).ranges({ data: "a" });
  const doc = b.getDocument();
  const r = await compile(doc, reg, fixtureLoaders());
  expect(r.issues).toEqual([]);
  set = r.set;
  view = resolveView(doc, set, "v", reg, supportedAll(reg)).view;
}, 60_000);

const annotation = () => set.annotations[view.data[0]] as NonNullable<SeqSet["annotations"][number]>;

describe("resolve's per-row paths allocate nothing", () => {
  it("cover: a row's codons marked through the map", () => {
    const path = view.paths[view.placements.path[0]];
    expect(path.span).toBe(3);
    const a = annotation();
    const ends: End[] = Array.from({ length: ROWS }, (_, i) => ({ system: path.from, slot: 0, side: 0, rows: [i] }));
    const mask = new Uint8Array(900);
    const step = (i: number) => cover(a, ends[i], path, mask);
    const r = heapScaling(N, step);
    expect(mask.indexOf(0)).toBe(-1);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    // The template's control (CF-05): the same step with a leak is caught.
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });

  it("the filter's test of a row, a string and a number field", () => {
    const a = annotation();
    const keep = matcher(["helix", 0.5, "1"]);
    let kept = 0;
    const r = heapScaling(N, (i) => {
      kept += +keep(a, set, i, "category") + +keep(a, set, i, "value");
    });
    expect(kept).toBeGreaterThan(0);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  // A continuous palette's color of a row is the viewer's per row (LY3 reads the number column): a double passed to
  // the palette's function is boxed, so it is not guarded here.
  it("a row's color through the layer's categorical palette", () => {
    const a = annotation();
    const f = layerPalette(set, [a], { kind: "categorical" }, "category", [0, 1]);
    let colored = 0;
    const r = heapScaling(N, (i) => {
      colored += f(paletteValue(a, set, i, "category")) === undefined ? 0 : 1;
    });
    expect(colored % N).toBe(0); // every row colored, in every run
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });
});
