/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A `pairs` arc's shape changes continuously as the view pans (seqspec:format#layer-types `pairs`): panned a px at a
// time, across the moments a partner enters and leaves the window and the widest arc enters and leaves the view, every
// arc's points move a px or so, never jump. The bug: an arc wider than the window snapped between forms (rising to the
// band's top past the edge, a line along the top once both ends were off) and the arcs' scale (the widest span in view)
// jumped as the widest arc came and went, so a short arc among long ones dropped or rose at once.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { case7, case21 } from "../../core/_test/format/cases.ts";
import type { PairsState } from "../src/prepare.ts";
import { Arcs } from "../src/renderer-helpers.ts";
import { LayerFrame, RecordingBackend } from "../src/testing/index.ts";
import { NO_SIZE } from "../src/themes.ts";
import { drawn } from "./draw-fixtures.ts";

/** Samples along each arc: at t = i / SAMPLES of its ends' span. */
const SAMPLES = 20;

/**
 * `doc`'s section at the zoom showing [start, stop) of `block`, scrolled to `from` px, panned right a px at a time
 * `steps` times (or to the axis's end): at each step, every pairs arc's
 * points in view (each sample's x and height, by layer, pair and sample) against the step before. Its largest moves.
 */
async function sweep(
  doc: SeqViewSpecDoc,
  view: string,
  section: string,
  zoom: [number, number, number],
  from: number,
  steps: number,
) {
  const t = await drawn(doc);
  t.ctx.views.setCurrent(view);
  await t.ctx.graph.update();
  for (const ref of t.ctx.viewports.sections.value) t.ctx.viewports.get(ref)?.setSize(1000, 400, 1);
  await t.ctx.graph.update();
  const vp = t.ctx.viewports.get(section);
  vp?.zoomTo(zoom[0], zoom[1], zoom[2], { durationMs: 0 });
  vp?.scrollTo(from, 0);
  await t.ctx.graph.update();
  const m = t.ctx.graph.assembled.value?.view;
  if (!m) throw new Error("no view");
  const entry = t.ctx.registry.get("renderers", "pairs");
  if (!entry) throw new Error("no pairs renderer");
  const frames: LayerFrame[] = [];
  for (let k = 0; k < m.prepared.length; k++) {
    const r = m.resolved;
    if (m.prepared[k].renderer !== "pairs") continue;
    if (m.set.strings.get(r.sections.ref[r.tracks.section[r.layers.track[k]]]) !== section) continue;
    frames.push(new LayerFrame(t.ctx, entry, section, k));
  }
  expect(frames.length).toBeGreaterThan(0);
  const rec = new RecordingBackend();
  let last = new Map<string, [number, number]>();
  let moveX = 0;
  let moveH = 0;
  let where = "";
  let compared = 0;
  let at = Number.NaN;
  for (let step = 0; step <= steps; step++) {
    if (step > 0) vp?.panBy(1, 0, "api");
    const now = new Map<string, [number, number]>();
    for (const f of frames) {
      rec.reset();
      f.draw(rec);
    }
    // At the axis's end the pan stops.
    if (frames[0].v.blockX[0] === at) break;
    at = frames[0].v.blockX[0];
    for (const f of frames) {
      const v = f.v;
      const s = m.prepared[f.k].state as PairsState;
      for (let k = 0; k < s.count; k++) {
        if (!Arcs.ends(v, s, k)) continue;
        const [xa, xb] = [s.xs[0], s.xs[1]];
        Arcs.apex(s, NO_SIZE);
        for (let i = 0; i <= SAMPLES; i++) {
          const x = xa + ((xb - xa) * i) / SAMPLES;
          if (x < v.left || x > v.right) continue;
          s.frame[10] = x;
          Arcs.heightAt(s);
          now.set(`${f.k}:${k}:${i}`, [x, s.frame[11]]);
        }
      }
    }
    for (const [key, [x, h]] of now) {
      const before = last.get(key);
      if (!before) continue;
      compared++;
      const ax = Math.abs(x - before[0]);
      const ah = Math.abs(h - before[1]);
      if (ax > moveX) moveX = ax;
      if (ah > moveH) {
        moveH = ah;
        where = `step ${step}, ${key}: ${before[1].toFixed(2)} → ${h.toFixed(2)}`;
      }
    }
    last = now;
  }
  return { moveX, moveH, where, compared };
}

/** A px a step: no sample moves more than 2 px across or up. */
async function smooth(p: Promise<Awaited<ReturnType<typeof sweep>>>) {
  const { moveX, moveH, where, compared } = await p;
  expect(compared).toBeGreaterThan(100);
  expect(moveX).toBeLessThanOrEqual(2);
  expect(moveH, where).toBeLessThanOrEqual(2);
}

describe("pairs: arcs move continuously as the view pans a px at a time", () => {
  it("case 21's sashimi on the gene's own strand, between zooms, the whole axis", async () => {
    await smooth(sweep(case21(), "gene2", "gene2.own", [0, 300, 700], 0, 3000));
  });
  it("case 21's sashimi on the gene's own strand at letter zoom, across exon 2's junctions", async () => {
    await smooth(sweep(case21(), "gene2", "gene2.own", [0, 734, 790], 1500, 2500));
  });
  it("case 21's sashimi on the region as written, between zooms, the whole axis", async () => {
    await smooth(sweep(case21(), "gene2", "gene2.reference", [0, 500, 900], 0, 3000));
  });
  it("case 7's contacts, zoomed in, from barnase into barstar", async () => {
    await smooth(sweep(case7(), "complex-one", "complex-one.main", [1, 30, 80], 0, 3000));
  });
});
