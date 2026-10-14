/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A viewer on the React UI for the browser heap guards (react.heap.ts), imported afresh per case: the default arrangement
// on the stub surface with the row controls on, a pointer moving inside one residue's cell (so a step is the UI's work
// for a pointer move and nothing else: the tooltip placed, the track's box asked), and a viewport whose version is bumped
// in place (the UI's work for a scroll, without the viewport's own).

import { act, Profiler } from "react";
import { createRoot } from "react-dom/client";
import "../../src/skin.css";
import { Inputs, PointerType } from "../../src/index.js";
import { SeqViewerLayout, SeqViewerProvider } from "../../src/react/index.js";
import { at, openedContext, stubSurface, WIDTH, xOf } from "../dom/fixtures.js";

(globalThis as { IS_REACT_ACT_ENVIRONMENT?: boolean }).IS_REACT_ACT_ENVIRONMENT = true;

export async function reactFixture(view: "chains" | "msa" = "chains") {
  const { ctx } = await openedContext();
  if (view === "msa") ctx.views.setCurrent("msa");
  const host = document.createElement("div");
  host.style.width = `${WIDTH}px`;
  document.body.appendChild(host);
  let renders = 0;
  const root = createRoot(host);
  const Control = () => <button type="button">3D</button>;
  await act(async () =>
    root.render(
      <Profiler id="ui" onRender={() => renders++}>
        <SeqViewerProvider context={ctx} ui={{ surface: stubSurface, components: { RowControl: Control } }}>
          <SeqViewerLayout />
        </SeqViewerProvider>
      </Profiler>,
    ),
  );
  await act(async () => new Promise((r) => setTimeout(r, 50)));
  const S = `${view}.main`;
  const p = Inputs.pointer();
  p.section = ctx.layout.sectionIndex(S);
  p.type = PointerType.move;
  const y = at(ctx, S, view === "msa" ? "t.rows" : "t.res").y;
  // Moves inside one residue's cell: the hovered element stays.
  const cell = xOf(ctx, S, 13) - xOf(ctx, S, 12);
  const x0 = xOf(ctx, S, 12) - cell / 2 + 0.5;
  const moveStep = (i: number) => {
    p.x = x0 + (i % 8) * ((cell - 1) / 8);
    p.y = y;
    ctx.input.pointer(p);
  };
  await act(async () => moveStep(0));
  // The viewport's version bumped in place, its scroll one of a few positions within a row (so the rows in view stay): what
  // the UI does for a scroll frame.
  const vp = ctx.viewports.get(S);
  const live = vp?.live as unknown as { current: { scrollY: number }; bump(): void };
  const scrollStep = (i: number, base = 0) => {
    live.current.scrollY = base + (i % 12) * 0.25;
    live.bump();
  };
  return {
    ctx,
    host,
    moveStep,
    scrollStep,
    renders: () => renders,
    dispose: async () => {
      await act(async () => root.unmount());
      host.remove();
      ctx.dispose();
    },
  };
}
