/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A header's hover is hover-rate (seqspec:interaction#hover, seqspec:viewer#performance): the hover alternating between
// two track headers and a group's header, each landing described by the `headerHover` behavior (its items counted, its
// runs written) and projected by the channels, allocates nothing after warm-up. Not vacuous: the described count
// alternates as the headers' items do, and the channels rebuilt the hover every step. The positive control: the same
// steps allocating an object each (`leaking`), warmed and measured alike, are caught. `toLoci` (item loci, allocating)
// is read only by a host, never here. A fresh module per case (`vi.resetModules`), as the draw guards are.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { expect, it, vi } from "vitest";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "./heap.ts";

const N = 50_000;
const S = "v.main";

function doc() {
  const b = createSeqViewSpec({ title: "Headers" });
  const a = b.data.sequence({ ref: "a", value: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ".repeat(4) });
  const ranges = (ref: string, n: number, step: number) => {
    const start = Array.from({ length: n }, (_, i) => step * i + 1);
    return b.data.annotation({
      ref,
      geometry: "ranges",
      target: a,
      value: { start, stop: start.map((s) => s + 1), label: start.map((s) => `${ref}${s}`) },
    });
  };
  const v = b.view({ ref: "v", label: "V" }).block({ sequence: a });
  v.track({ ref: "v.forty", label: "Forty" }).ranges({ data: ranges("d40", 40, 3) });
  v.track({ ref: "v.ten", label: "Ten" }).ranges({ data: ranges("d10", 10, 7) });
  v.group({ ref: "v.g", label: "Group", collapsed: false } as never)
    .track({ ref: "v.g.five", label: "Five" })
    .ranges({ data: ranges("d5", 5, 11) });
  return b.getDocument();
}

/** The headers' document opened afresh, and a step: the hover on a header by turns, described and projected. */
async function headers() {
  vi.resetModules();
  const { opened } = await import("./layout-fixtures.ts");
  const { ChannelMarks, EntryKind, Hit } = await import("../src/index.ts");
  const t = await opened(doc(), { size: { [S]: [1000, 500] } });
  t.ctx.behaviors.get("hoverLabel")?.remove();
  const l = t.layout(S);
  const at = (y: number) => {
    const h = Hit.create();
    t.ctx.layout.hitTest(S, 10.5, y, h);
    return h;
  };
  let groupY = -1;
  for (let e = 0; e < l.entryCount; e++)
    if (l.entryKind[e] === EntryKind.group) groupY = l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2;
  const hits = [at(t.yOf(S, "v.forty")), at(t.yOf(S, "v.ten")), at(groupY)];
  const ch = new ChannelMarks(t.ctx);
  const header = t.ctx.hover.live.current.header;
  const counts = [0, 0, 0];
  const out = { t, ch, counts, sink: 0, step: (_: number) => {} };
  out.step = (i: number) => {
    const k = i % 3;
    t.ctx.hover.setFromHit(hits[k], "pointer");
    ch.sync();
    counts[k] = header.items;
    out.sink += ch.marks(0)?.hover.count ?? 0;
  };
  return out;
}
/** The guard's measurement: 6,000 steps warm, `warmed` called, then `N` and `4N` measured. */
function measure(step: (i: number) => void, warmed = () => {}) {
  for (let i = 0; i < 6_000; i++) step(i);
  warmed();
  return heapScaling(N, step);
}

it("headers alternating, each described and projected: flat; the count alternates, the channels rebuild", async () => {
  const h = await headers();
  let before = 0;
  const r = measure(h.step, () => {
    before = h.ch.rebuilds.hover;
  });
  expect(isFlat(r), describeScaling(r, N)).toBe(true);
  expect(h.counts).toEqual([40, 10, 5]);
  expect(h.ch.rebuilds.hover - before).toBeGreaterThanOrEqual(5 * N);
  expect(h.sink).toBeGreaterThan(0);
  h.t.ctx.dispose();
}, 120_000);

it("the same steps allocating an object each are caught (positive control)", async () => {
  const h = await headers();
  expect(controlCaught(() => measure(leaking(h.step)))).toBe(true);
  h.t.ctx.dispose();
}, 120_000);
