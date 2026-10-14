/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guard for LY3b's hits (seqspec:viewer#performance): a hit test per pointer move over the kernel entries' own
// hits (letters `of` another polymer, ranges in lanes, pairs' arcs and ends over hundreds of pairs, a matrix's state
// bands), on `heapScaling` (a step one move). The positive control: the same moves allocating an object each
// (`leaking`), through the same measurement, are caught.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { expect, it } from "vitest";
import { case14 } from "../../core/_test/format/cases.ts";
import { Hit, HoverKind, KernelRenderers } from "../src/index.ts";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "./heap.ts";
import { opened } from "./layout-fixtures.ts";
import { allTypes, rangesDoc } from "./prepare-fixtures.ts";

const N = 50_000;
let sink = 0;

function matrixDoc() {
  const b = createSeqViewSpec({ title: "Matrix" });
  const a = b.data.sequence({ ref: "a", label: "A", value: "MKTAYIAKQR".repeat(10) });
  const index = Array.from({ length: 100 }, (_, i) => i);
  const value = { index, value: index.map((i) => [i % 3, (i * 7) % 5, 1]) };
  const m = b.data.annotation({ ref: "m", geometry: "matrix", target: a, states: ["x", "y", "z"], value } as never);
  b.view({ ref: "v", label: "V" }).block({ sequence: a }).track({ label: "M" }).matrix({ data: m });
  return b.getDocument();
}

it("hit tests per move over the kernel entries' hits allocate nothing", async () => {
  const docs = [
    allTypes({ length: 2000 }),
    case14(),
    rangesDoc([
      [10, 30],
      [15, 20],
      [12, 40],
    ]),
    matrixDoc(),
  ];
  const opens = await Promise.all(docs.map((d) => opened(d, { renderers: [...KernelRenderers] })));
  const points: [(typeof opens)[number], number, number, number][] = [];
  for (const t of opens) {
    const s = t.ctx.viewports.sections.value[0];
    t.vp(s).setView(4, 0, "api");
    const l = t.layout(s);
    const i = t.ctx.layout.sectionIndex(s);
    for (let y = l.navHeight; y < Math.min(l.height, l.scrolledTop + 400); y += 3)
      for (let x = l.dataLeft + 1; x < l.dataRight; x += 37) points.push([t, i, x, y]);
  }
  const hit = Hit.create();
  let items = 0;
  let states = 0;
  for (const p of points) {
    p[0].ctx.layout.hitTest(p[1], p[2], p[3], hit);
    if (hit.kind === HoverKind.item) items++;
    if (hit.state >= 0) states++;
  }
  // Not vacuous: items and matrix states are hit.
  expect([items > 20, states > 20]).toEqual([true, true]);
  const step = (k: number) => {
    const p = points[k % points.length];
    p[0].ctx.layout.hitTest(p[1], p[2], p[3], hit);
    sink = (sink + hit.kind + hit.item + hit.state) & 0xffff;
  };
  const r = heapScaling(N, step);
  expect(isFlat(r), describeScaling(r, N)).toBe(true);
  expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
});
