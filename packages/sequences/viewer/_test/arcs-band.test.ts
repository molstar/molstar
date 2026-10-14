/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Every `pairs` arc stays within its track's band at every zoom (seqspec:format#layer-types `pairs`): fitted, between
// and at letter zoom, case 21's sashimi on the gene's own strand and on the region as written (arcs hanging below and
// rising above their ends) and case 7's contacts. The bug: an arc wider than the window had its apex past the band, so
// zoomed in its visible piece was a steep, near-straight line out of the band into the next section.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { case7, case21 } from "../../core/_test/format/cases.ts";
import { type PathBatch, PathMode, PathVerb } from "../src/draw/backend.ts";
import { LayerFrame, RecordingBackend } from "../src/testing/index.ts";
import { drawn } from "./draw-fixtures.ts";

/** Points along each subpath of a stroked batch: its moves and lines, and 11 along each quad. */
function points(b: PathBatch): [number, number][] {
  const out: [number, number][] = [];
  for (let i = 0; i < b.n; i++) {
    if (b.mode[i] !== PathMode.stroke) continue;
    const v1 = i + 1 < b.n ? b.start[i + 1] : b.nv;
    let o = b.cstart[i];
    let x = 0;
    let y = 0;
    for (let j = b.start[i]; j < v1; j++) {
      const verb = b.verbs[j];
      if (verb === PathVerb.quad) {
        const [cx, cy, ex, ey] = [b.coords[o], b.coords[o + 1], b.coords[o + 2], b.coords[o + 3]];
        for (let q = 0; q <= 10; q++) {
          const t = q / 10;
          const [a, m, c] = [(1 - t) * (1 - t), 2 * t * (1 - t), t * t];
          out.push([a * x + m * cx + c * ex, a * y + m * cy + c * ey]);
        }
        [x, y] = [ex, ey];
        o += 4;
      } else {
        [x, y] = [b.coords[o], b.coords[o + 1]];
        out.push([x, y]);
        o += 2;
      }
    }
  }
  return out;
}

/** Each `pairs` layer of `section`, drawn alone at the section's view: its band and its arcs' points. */
async function arcsAt(doc: SeqViewSpecDoc, view: string, section: string, zoom: [number, number, number] | null) {
  const t = await drawn(doc);
  t.ctx.views.setCurrent(view);
  await t.ctx.graph.update();
  for (const ref of t.ctx.viewports.sections.value) t.ctx.viewports.get(ref)?.setSize(1000, 400, 1);
  await t.ctx.graph.update();
  if (zoom) t.ctx.viewports.get(section)?.zoomTo(zoom[0], zoom[1], zoom[2], { durationMs: 0 });
  await t.ctx.graph.update();
  const m = t.ctx.graph.assembled.value?.view;
  if (!m) throw new Error("no view");
  const out: { top: number; bottom: number; pts: [number, number][] }[] = [];
  for (let k = 0; k < m.prepared.length; k++) {
    const r = m.resolved;
    if (m.prepared[k].renderer !== "pairs") continue;
    if (m.set.strings.get(r.sections.ref[r.tracks.section[r.layers.track[k]]]) !== section) continue;
    const entry = t.ctx.registry.get("renderers", "pairs");
    if (!entry) throw new Error("no pairs renderer");
    const f = new LayerFrame(t.ctx, entry, section, k);
    const rec = new RecordingBackend();
    f.draw(rec);
    out.push({ top: f.v.top, bottom: f.v.top + f.v.height, pts: points(rec.pathItems) });
  }
  return out;
}

const ZOOMS: [string, [number, number, number] | null][] = [
  ["fitted", null],
  ["mid", [0, 300, 700]],
  ["letters 735–790", [0, 734, 790]],
  ["letters 600–650", [0, 599, 650]],
];

describe("pairs: every arc within its track's band", () => {
  for (const section of ["gene2.own", "gene2.reference"])
    for (const [name, zoom] of ZOOMS)
      it(`case 21's sashimi, ${section}, ${name}`, async () => {
        const layers = await arcsAt(case21(), "gene2", section, zoom);
        expect(layers.length).toBe(2);
        let n = 0;
        let out = 0;
        for (const { top, bottom, pts } of layers) {
          for (const [, y] of pts) if (y < top || y > bottom) out++;
          n += pts.length;
        }
        expect(n).toBeGreaterThan(0);
        expect(out).toBe(0);
      });

  for (const [name, zoom] of [
    ["fitted", null],
    ["barnase 20–70", [0, 20, 70]],
    ["barstar 30–80", [1, 30, 80]],
    ["barstar 36–60", [1, 36, 60]],
  ] as const)
    it(`case 7's contacts, ${name}`, async () => {
      const layers = await arcsAt(case7(), "complex-one", "complex-one.main", zoom as [number, number, number] | null);
      expect(layers.length).toBe(1);
      const [{ top, bottom, pts }] = layers;
      expect(pts.length).toBeGreaterThan(0);
      expect(pts.filter(([, y]) => y < top || y > bottom)).toEqual([]);
    });
});
