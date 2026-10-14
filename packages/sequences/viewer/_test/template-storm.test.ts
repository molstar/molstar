/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1b-b's templated storm: 10,000 dynamic nodes held, 1,000 of them under a template's source (1,000 instance tracks),
// 9,000 under another source. Then one node a frame: under the template's source (an arrival: its template and the view
// rerun, every other instance kept), under the other (only `#assemble.stream`), and evictions of instances. Reports the
// times and, apart, what an update spends in the set's `withChanges` and in `templateNodes` at 10,000 held; asserts
// what reruns and loose times. The slow tier's; the fast tier's storm of 50 instances over 500 checks what reruns.

import { indexDoc, templateNodes } from "@molstar/seq-core/format";
import { SeqSets } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { openTemplates, P, range, templateDoc } from "./template-fixtures.ts";

const median = (xs: number[]) => xs.slice().sort((a, b) => a - b)[xs.length >> 1];
const ms = (x: number) => Math.round(x * 1000) / 1000;
const at = (k: number) => 1 + (k % (P.length - 8));

/**
 * `instances` nodes under the template's source and the rest of `nodes` under another, held; then `frames` of one node
 * a frame arriving under each and leaving. What reran and the instances at each step are checked; the times returned.
 */
async function storm(instances: number, nodes: number, frames: number) {
  const t = await openTemplates({});
  const m = t.manager;
  let t0 = performance.now();
  for (let k = 0; k < instances; k++) m.nodes.add("genes", [range(`g${k}`, at(k))]);
  for (let k = 0; k < nodes - instances; k++) m.nodes.add("notes", [range(`n${k}`, at(k))]);
  const hold = performance.now() - t0;
  t0 = performance.now();
  await t.frame();
  const once = performance.now() - t0;
  expect(t.instances()).toHaveLength(instances);
  expect(m.nodes.entries.size).toBe(nodes);

  const timed = async (fn: (k: number) => void) => {
    const each: number[] = [];
    const reran = new Set<string>();
    for (let k = 0; k < frames; k++) {
      fn(k);
      const s = performance.now();
      await t.frame();
      each.push(performance.now() - s);
      for (const r of t.ctx.graph.ran) reran.add(r);
    }
    return { each, reran: [...reran].sort() };
  };
  const arrive = await timed((k) => m.nodes.add("genes", [range(`new${k}`, at(k))]));
  expect(arrive.reran).toEqual(["#assemble.stream", "#template:tpl", "#view:v"]);
  expect(t.instances()).toHaveLength(instances + frames);
  const other = await timed((k) => m.nodes.add("notes", [range(`more${k}`, at(k))]));
  expect(other.reran).toEqual(["#assemble.stream"]);
  const leave = await timed((k) => m.nodes.release([`genes/new${k}`]));
  expect(leave.reran).toEqual(["#assemble.stream", "#template:tpl", "#view:v"]);
  expect(t.instances()).toHaveLength(instances);
  return { t, hold, once, arrive, other, leave };
}

describe("a templated storm", () => {
  it("50 instances over 500 held nodes: an arrival reruns its template and the view, another source's nothing", async () => {
    const { t } = await storm(50, 500, 10);
    t.ctx.dispose();
  });

  it("1,000 instances over 10,000 held nodes: an arrival reruns its template and the view, another source's nothing", {
    tags: ["slow"],
    timeout: 300_000,
  }, async () => {
    const { t, hold, once, arrive, other, leave } = await storm(1_000, 10_000, 200);

    // Apart: one annotation upserted into the set holding 10,000, and a template's selection over it.
    const set = t.ctx.graph.assembled.value?.set;
    if (!set) throw new Error("no set");
    const one = { key: "genes/probe", geometry: "ranges" as const, target: "p", columns: { start: [1], stop: [4] } };
    const w: number[] = [];
    for (let k = 0; k < 50; k++) {
      const s = performance.now();
      SeqSets.withChanges(set, { upsert: { annotations: [one] } });
      w.push(performance.now() - s);
    }
    const index = indexDoc(templateDoc(), t.ctx.registry.svs);
    const sel: number[] = [];
    for (let k = 0; k < 50; k++) {
      const s = performance.now();
      templateNodes("tpl", index, set);
      sel.push(performance.now() - s);
    }
    const line = (x: { each: number[] }) => `median ${ms(median(x.each))} ms, max ${ms(Math.max(...x.each))} ms`;
    process.stdout.write(
      `ST1b-b templated storm: hold 10,000 ${ms(hold)} ms; first update (1,000 instances) ${ms(once)} ms; ` +
        `an arrival under the template ${line(arrive)}; under another source ${line(other)}; an eviction ${line(leave)}; ` +
        `withChanges at 10,000 median ${ms(median(w))} ms; templateNodes median ${ms(median(sel))} ms\n`,
    );
    expect(median(arrive.each)).toBeLessThan(50);
    expect(median(other.each)).toBeLessThan(25);
    t.ctx.dispose();
  });
});
