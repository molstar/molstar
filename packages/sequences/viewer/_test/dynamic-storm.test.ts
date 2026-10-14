/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1b-a's patch storm: 1,000 node patches over the streamed genome with the dynamic tracks, held straight into the
// stream manager's nodes (one node a patch, as a region's would bring it), then applied: all in one frame (one graph
// update), and one a frame (1,000 updates at stream rate, the case the assembly must not rebuild the view for), then
// the same 1,000 again (no-ops). Reports the time a patch costs to hold, an update's time and what it reran; asserts
// what reruns, and loose times. The slow tier's; the fast tier's storm of 50 checks what reruns.

import { describe, expect, it } from "vitest";
import { GRCH38 } from "../../core/_test/format/stream-provider.ts";
import type { Streams } from "../src/streams.ts";
import { dynamicDoc, GENES, openGenome, peakNode, type StreamedGenome, withPeaks } from "./stream-fixtures.ts";

const node = (k: number) => peakNode(GRCH38[k % 3][0], 4096 * (1000 + k));
const median = (xs: number[]) => xs.slice().sort((a, b) => a - b)[xs.length >> 1];

async function frame(t: StreamedGenome) {
  t.scheduler.advance(0);
  await t.ctx.graph.update();
}

/** `n` node patches held then applied in one frame, `n` more one a frame, then the first `n` again; the times. */
async function storm(n: number) {
  const t = await openGenome({
    chromosomes: GRCH38.slice(0, 3),
    doc: dynamicDoc(),
    nodes: [GENES],
    answer: withPeaks(),
  });
  const m = t.manager as Streams;
  const set0 = t.ctx.graph.assembled.value?.set;
  const view0 = t.ctx.graph.assembled.value?.view;

  // All at once: n patches held, then one frame applies them.
  let t0 = performance.now();
  for (let k = 0; k < n; k++) m.nodes.add("phylop_src", [node(k)]);
  const hold = (performance.now() - t0) / n;
  t0 = performance.now();
  await frame(t);
  const once = performance.now() - t0;
  const set1 = t.ctx.graph.assembled.value?.set;
  expect(m.nodes.entries.size).toBe(n);
  expect(set1?.lineage).toBe(set0?.lineage);

  // At stream rate: one patch a frame, n frames over the n held.
  const each: number[] = [];
  const reran = new Set<string>();
  for (let k = 0; k < n; k++) {
    m.nodes.add("phylop_src", [node(n + k)]);
    t0 = performance.now();
    await frame(t);
    each.push(performance.now() - t0);
    for (const r of t.ctx.graph.ran) reran.add(r);
  }
  // A node no track names changes the set only: no section, layer or view reads it (ST1b-b).
  expect([...reran].sort()).toEqual(["#assemble.stream"]);
  expect(t.ctx.graph.assembled.value?.view?.resolved.layers.ref).toEqual(view0?.resolved.layers.ref);

  // The same n again: no-ops, the set and the version kept.
  const version = m.nodes.version;
  const set2 = t.ctx.graph.assembled.value?.set;
  t0 = performance.now();
  for (let k = 0; k < n; k++) m.nodes.add("phylop_src", [node(k)]);
  const again = (performance.now() - t0) / n;
  await frame(t);
  expect([m.nodes.version, t.ctx.graph.assembled.value?.set]).toEqual([version, set2]);
  t.ctx.dispose();
  return { hold, once, each, again };
}

describe("a storm of dynamic nodes", () => {
  it("50 node patches: applied in one update a frame that reruns no section, layer or view; again, no-ops", async () => {
    await storm(50);
  });

  it("1,000 node patches: held in microseconds each, applied in one update a frame that reruns no section, layer or view", {
    tags: ["slow"],
    timeout: 120_000,
  }, async () => {
    const { hold, once, each, again } = await storm(1_000);
    const ms = (x: number) => Math.round(x * 1000) / 1000;
    process.stdout.write(
      `ST1b-a storm: hold ${ms(hold * 1000)} µs a patch (again ${ms(again * 1000)} µs); 1,000 in one update ${ms(once)} ms; ` +
        `one a frame over 1,000 to 2,000 held: median ${ms(median(each))} ms, max ${ms(Math.max(...each))} ms\n`,
    );
    expect(median(each)).toBeLessThan(16);
  });
});
