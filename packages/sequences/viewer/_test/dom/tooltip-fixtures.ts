/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A viewer with its tooltip shown, for the tooltip's browser heap guards (tooltip.heap.ts), imported afresh per case:
// the default UI on the stub surface, a pointer moving inside one residue's cell so that moves change only the
// pointer, and a keyboard cursor whose viewport moves.

import * as D from "../../src/dom/index.ts";
import { Inputs, PointerType } from "../../src/index.ts";
import "../../src/skin.css";
import { at, openedContext, stubSurface, WIDTH, xOf } from "./fixtures.ts";

export async function tooltipFixture() {
  const { ctx } = await openedContext();
  const host = document.createElement("div");
  host.style.width = `${WIDTH}px`;
  document.body.appendChild(host);
  const mounted = D.mountSeqViewer(host, ctx, { surface: stubSurface });
  const S = "chains.main";
  const y = at(ctx, S, "t.res").y;
  const p = Inputs.pointer();
  p.section = ctx.layout.sectionIndex(S);
  p.type = PointerType.move;
  // Moves inside one residue's cell: the hovered element stays, so a step is the pointer's move and the tooltip's placement.
  const cell = xOf(ctx, S, 13) - xOf(ctx, S, 12);
  const x0 = xOf(ctx, S, 12) - cell / 2 + 0.5;
  const moveStep = (i: number) => {
    p.x = x0 + (i % 8) * ((cell - 1) / 8);
    p.y = y;
    ctx.input.pointer(p);
  };
  moveStep(0);
  let hovered = 0;
  ctx.hover.live.version.subscribe(() => hovered++);
  const tip = host.querySelector(".seqv-tooltip") as HTMLElement;
  if (!tip) throw new Error("no tooltip");
  const cursorStep = () => {
    ctx.input.focus(p.section, true);
    // Keep the cursor's element fixed and bump the viewport's version instead: the tooltip is placed again from it.
    const vp = ctx.viewports.get(S);
    return (i: number) => vp?.panBy((i & 1) === 0 ? 1 : -1, 0, "pointer");
  };
  return {
    ctx,
    tip,
    moveStep,
    cursorStep,
    hovers: () => hovered - 1,
    dispose() {
      mounted.dispose();
      host.remove();
      ctx.dispose();
    },
  };
}
