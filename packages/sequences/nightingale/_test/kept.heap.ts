/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guard for LS1's per-hover path in `mount`: a section's `highlights` (`AxisMeta.highlights`) kept under a host's
// highlight (`keepLit`) on every change, after the host's loci are projected, on RB0b's stepped `heapScaling`. The
// positive control is the template's: the same step with a leak, at the guard's own length (NG-11: the strings it kept on
// a list could cancel in the shorter run).

import { projectionOf } from "@molstar/seq-core/format";
import { Extents, type Loci, OrderedSet, type UnitId } from "@molstar/seq-core/model";
import { beforeAll, describe, expect, it } from "vitest";
import { gene, open } from "../../core/_test/format/project-fixtures.ts";
import { controlCaught, describeScaling, heapScaling, isFlat, leak } from "../../core/_test/model/heap.ts";
import { AxisRuns, AxisSink, axesOf, highlightRuns, keepLit, prepareLoci } from "../src/axis.ts";
import type { AxisMeta } from "../src/index.ts";

const N = 50_000;
let g: Awaited<ReturnType<typeof open>>;
let meta: AxisMeta;
let loci: Loci[];

beforeAll(async () => {
  g = await open(gene({ chromosome: true }));
  const kept = AxisRuns.create();
  for (let k = 0; k < 40; k++) AxisRuns.add(kept, 7 * k, 7 * k + 3);
  meta = { ...axesOf(g.view, g.set, 0)[0], highlights: AxisRuns.finish(kept) };
  const unit = g.system("gene/unit").id as UnitId;
  loci = Array.from({ length: 120 }, (_, i) => ({
    kind: "seq-loci",
    set: g.set,
    elements: [{ unit, indices: OrderedSet.ofSingleton(i) as never }],
  }));
}, 60_000);

describe("a section's highlights kept under a host's", () => {
  const run = (n: number, also: (out: AxisRuns) => void) => {
    const [p, ext, sink, out] = [projectionOf(g.view, g.set), Extents.create(8), new AxisSink(), AxisRuns.create()];
    let lit = 0;
    const r = heapScaling(n, (i) => {
      if (prepareLoci(p, ext, loci[i % loci.length], "element")) highlightRuns(p, ext, meta, sink, out);
      else out.count = 0;
      lit += keepLit(out, meta.highlights).count;
      also(out);
    });
    expect(lit).toBeGreaterThan(40 * n);
    return r;
  };

  it("allocate nothing per change", () => {
    const r = run(N, () => {});
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("the same step with a leak is caught (the positive control)", () => {
    expect(controlCaught(() => run(N, (out) => leak(out.count)))).toBe(true);
  });
});
