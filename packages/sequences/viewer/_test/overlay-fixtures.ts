/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for the composer's overlay (CS1b): documents with stacked blocks, wrapped lines and an alignment's rows
// (gaps, data columns), loci on them, a recording's rects by pass, and hover frames for the heap guards
//, which import this module afresh per case.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { type Loci, OrderedSet, type SeqSet } from "@molstar/seq-core/model";
import { Dirty } from "../src/draw/index.js";
import { Hit, Inputs, PointerType } from "../src/index.js";
import { Op, type RecordingBackend } from "../src/testing/index.js";
import { COPY_BACK, type ComposedOptions, composed, recording } from "./composer-fixtures.js";
import { directed } from "./conformance-fixtures.js";
import { family, rowLoci, seqLoci } from "./overlay-docs.js";
import { scaleDoc } from "./prepare-fixtures.js";

export * from "./overlay-docs.js";

/** A recording's rects in pass `pass` as [x, y, w, h, color, alpha, stroke]. */
export function passRects(r: RecordingBackend, pass: number): number[][] {
  const b = r.rectItems;
  const out: number[][] = [];
  let at = -1;
  for (let k = 0; k < r.opCount; k++) {
    const code = r.ops[4 * k];
    if (code === Op.beginPass) at = r.args[r.ops[4 * k + 3]];
    if (code !== Op.rects || at !== pass) continue;
    for (let i = r.ops[4 * k + 1]; i < r.ops[4 * k + 1] + r.ops[4 * k + 2]; i++)
      out.push([b.x[i], b.y[i], b.w[i], b.h[i], b.color[i], b.alpha[i], b.stroke[i]]);
  }
  return out;
}

/** A context showing `doc` with the kernel renderers and its section's composer (composer-fixtures' `composed`). */
export async function lit(doc: SeqViewSpecDoc, o: ComposedOptions = {}) {
  const t = await composed(doc, o);
  const set = t.ctx.graph.assembled.value?.set as SeqSet;
  const hit = Hit.create();
  /** The hover on what is at (x, y), as the pointer's. */
  const hover = (x: number, y: number) => {
    if (t.ctx.layout.hitTest(t.section, x, y, hit)) t.ctx.hover.setFromHit(hit, "pointer");
  };
  const overlay = (into = recording()) => t.frame(into, Dirty.overlay);
  return { ...t, set, hover, overlay };
}

/**
 * Hover frames: the pointer moving over `scale` (2 × 200,000 positions zoomed out, a selection and a highlight on them)
 * or `rows` (an alignment's rows, a row's residues selected, cells highlighted), an overlay frame per move; `direction`
 * -1: every piece reversed (AX1).
 */
export async function hoverFrames(kind: "scale" | "rows", direction: 1 | -1 = 1) {
  const rows = kind === "rows";
  const doc = rows ? family() : scaleDoc();
  const t = await lit(direction < 0 ? directed(doc, -1) : doc, { width: 900, height: 420, dpr: 2 });
  const be = recording(COPY_BACK);
  // The overlay alone: the hover label's strings are one allocation per hovered element by design (IN2's guard).
  t.ctx.behaviors.get("hoverLabel")?.remove();
  if (rows) {
    t.vp?.setView(40, 0, "api");
    t.ctx.selection.set(rowLoci(t.set, 1, 0, 6), "api");
    t.ctx.highlight.set("host", {
      kind: "aligned-loci",
      set: t.set,
      alignment: 0,
      elements: [{ row: 2, columns: OrderedSet.ofRange(1, 4) }],
    } as unknown as Loci);
  } else {
    t.vp?.setView(0.0019, 0, "api");
    const units = Array.from(t.set.units.key, (k) => t.set.strings.get(k).replace("/unit", ""));
    t.ctx.selection.set(seqLoci(t.set, units[0], 1000, 90_000), "api");
    t.ctx.highlight.set("host", seqLoci(t.set, units[units.length - 1], 50_000, 150_000));
  }
  t.frame(be, Dirty.all);
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(t.section);
  p.type = PointerType.move;
  const l = t.ctx.layout.get(t.section)?.current;
  if (!l) throw new Error("no layout");
  const r = l.rows;
  const top = rows ? l.scrolledTop + l.entryTop[r.entry] + r.height / 2 : l.scrolledTop + 4;
  const frame = (i: number) => {
    p.x = rows ? l.dataLeft + 2 + (i % 50) * 4.5 : 200 + (i % 50) * 13;
    p.y = top + (i % (rows ? 4 : 9)) * (rows ? r.height : 30);
    t.ctx.input.pointer(p);
    t.frame(be, Dirty.overlay);
  };
  const channels = (t.composer as unknown as { channels: { rebuilds: { hover: number } } }).channels;
  return { frame, hovers: () => channels.rebuilds.hover };
}
