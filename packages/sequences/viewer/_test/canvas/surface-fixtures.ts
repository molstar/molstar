/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Frames through a mounted surface for the browser heap guards (surface.heap.ts), imported afresh per case: a scale
// document over two blocks with a selection and a highlight, a page-like scheduler holding a few frames (the manual one
// makes a Map per `advance`), and a step that sends events and runs the frame they asked for.

import { SectionSurface } from "../../src/canvas/index.ts";
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  Inputs,
  PointerType,
  type Scheduler,
  type SeqViewerContext,
} from "../../src/index.ts";
import { fixtureMeasurer } from "../../src/testing/index.ts";
import { seqLoci } from "../overlay-docs.ts";
import { scaleDoc } from "../prepare-fixtures.ts";

/**
 * Up to four outstanding frames, run in order by `flush` (the surface's and the dispatcher's re-hit, `InputDispatcher.moved`,
 * as a page's rAF queues both); timers never fire (a gesture's commit stays pending, as mid-gesture). No allocation.
 */
class OneFrame implements Scheduler {
  private readonly cbs: (((t: number) => void) | null)[] = [null, null, null, null];
  private readonly running: (((t: number) => void) | null)[] = [null, null, null, null];
  requested = 0;
  now() {
    return 0;
  }
  requestFrame(cb: (t: number) => void) {
    this.requested++;
    for (let k = 0; k < this.cbs.length; k++)
      if (this.cbs[k] === null) {
        this.cbs[k] = cb;
        return k + 1;
      }
    throw new Error("more than four frames outstanding");
  }
  cancelFrame(id: number) {
    if (id >= 1 && id <= this.cbs.length) this.cbs[id - 1] = null;
  }
  setTimer() {
    return 99;
  }
  clearTimer() {}
  flush() {
    for (let k = 0; k < this.cbs.length; k++) {
      this.running[k] = this.cbs[k];
      this.cbs[k] = null;
    }
    for (let k = 0; k < this.running.length; k++) {
      const cb = this.running[k];
      this.running[k] = null;
      cb?.(0);
    }
  }
}

export async function surfaceFrames(kind: "hover" | "pan" | "both") {
  const scheduler = new OneFrame();
  const ctx: SeqViewerContext = createSeqViewer(DefaultSeqViewerSpec(), { scheduler, text: fixtureMeasurer() });
  const r = await ctx.setDocument(scaleDoc());
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  await ctx.graph.update();
  const host = document.createElement("div");
  host.style.cssText = "position:absolute;left:0;top:0;width:900px;height:420px";
  document.body.appendChild(host);
  const section = ctx.viewports.sections.value[0];
  const surface = new SectionSurface(ctx, section, host, { dpr: 2 });
  // The surface's paths alone: the hover label's strings are one allocation per hovered element by design (IN2).
  ctx.behaviors.get("hoverLabel")?.remove();
  const vp = ctx.viewports.get(section);
  vp?.setView(kind === "pan" ? 0.4 : 0.0019, 0, "api");
  const set = ctx.graph.assembled.value?.set;
  if (!set) throw new Error("no set");
  ctx.selection.set(seqLoci(set, "chr0", 1000, 90_000), "api");
  ctx.highlight.set("host", seqLoci(set, "chr1", 50_000, 150_000));
  if (!ctx.selection.current.value.seq) throw new Error("the selection is not set");
  surface.renderNow();
  const p = Inputs.pointer();
  p.section = ctx.layout.sectionIndex(section);
  p.type = PointerType.move;
  const l = ctx.layout.get(section)?.current;
  if (!l) throw new Error("no layout");
  const top = l.scrolledTop + 4;
  const step = (i: number) => {
    if (kind !== "pan") {
      p.x = 200 + (i % 50) * 13;
      p.y = top + (i % 9) * 30;
      ctx.input.pointer(p);
    }
    if (kind !== "hover") vp?.panBy((i & 1) === 0 ? 9 : -9, 0, "pointer");
    scheduler.flush();
  };
  const dispose = () => {
    surface.dispose();
    ctx.dispose();
    host.remove();
  };
  return { step, surface, scheduler, dispose };
}
