/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The link's 3D-hover path (seqspec:libraries#the-mol-adapter, Rates): a Mol* loci to runs in the member systems and on
// to Seq* units, into reused `Extents`, allocates nothing, on RB0b's stepped guard. The positive control: the same step
// allocating an object each (`leaking`), measured alike, is caught. The Seq*-to-3D hover makes a Mol* loci per hover
// change by design (the link's `react`), so it has no guard.

import { compile } from "@molstar/seq-core/format";
import { CoordinateSystemKind, Extents, KeyKind } from "@molstar/seq-core/model";
import { Translator } from "@molstar/seq-core/translate";
import { describe, expect, it } from "vitest";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "../../core/_test/model/heap.ts";
import { readStructure, StructureLink, structureOf } from "../src/index.ts";
import { bytesOf, loaders, p53Doc, registries } from "./fixtures.ts";

const N = 50_000;
const { set } = await compile(p53Doc([{ ref: "1tup", pdb: "1TUP", entity: "3" }]), registries(), loaders());
const structure = await structureOf(await readStructure(bytesOf("1tup.bcif"), "bcif"), 0, null);
const link = StructureLink.attach(new Translator(), "3d", structure, set, "1tup");
const t = { kind: 0, id: 0, row: -1 };
set.byName("1tup/C", t);
expect(t.kind).toBe(KeyKind.unit);
const at = (i: number) => {
  const e = Extents.create(1);
  Extents.add(e, { kind: CoordinateSystemKind.unit, id: t.id, row: -1 }, i, i + 1);
  return link.toLoci(e);
};
// Hovered residues of chain C, made once (Mol* makes its hover loci; the link only reads them).
const hovers = Array.from({ length: 64 }, (_, i) => at(i * 3));
const members = Extents.create(8);
const seq = Extents.create(8);
let sink = 0;
/** A Mol* hover to member runs and Seq* runs. */
const hover = (i: number) => {
  members.count = 0;
  seq.count = 0;
  link.fromLoci(hovers[i & 63], members);
  link.toSeq(members, seq);
  sink = (sink + seq.start[0]) | 0;
};

describe("link heap guards", () => {
  it("a Mol* hover to member runs and Seq* runs", () => {
    const r = heapScaling(N, hover);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(sink).not.toBe(0.5);
  });

  it("the same step allocating an object each is caught (positive control)", () => {
    expect(controlCaught(() => heapScaling(N, leaking(hover)))).toBe(true);
  });
});
