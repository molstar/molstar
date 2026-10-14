/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Mounted sections for the browser tests (CS2's ports of the old stack's canvas fixtures): a context on the manual clock
// showing a document, its first section mounted in a host of a given size, the pixels of any canvas read back, and the
// comparison of two captures (premultiplied, the way the canvas stores them).

import type { Loaders, SeqViewSpecDoc } from "@molstar/seq-core/format";
import { mountSection, type SectionSurface } from "../../src/canvas/index.ts";
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  KernelRenderers,
  type LayerRenderer,
  manualScheduler,
  type RenderProps,
  type SeqViewerContext,
  type SeqViewerRegistryEntry,
} from "../../src/index.ts";
import { fixtureMeasurer } from "../../src/testing/index.ts";

export interface Opened {
  ctx: SeqViewerContext;
  surface: SectionSurface;
  host: HTMLDivElement;
  section: string;
  advance(ms?: number): void;
}

/** A host's tokens: fonts the draws are the same on any machine, the colors skin.css's. */
const HOST = '--seq-font-ui:Arial;--seq-font-mono:"Courier New"';

/** A context on the manual clock showing `doc`, its first section mounted in a `width` × `height` host. */
export async function opened(
  doc: SeqViewSpecDoc,
  o: {
    width?: number;
    height?: number;
    dpr?: number;
    render?: Partial<RenderProps>;
    registry?: SeqViewerRegistryEntry[];
    /** Frames and timers never fire (heap guards drive `renderNow` themselves). */
    still?: boolean;
    /** The label providers (default true): their strings are one allocation per hovered element. */
    labels?: boolean;
    loaders?: Loaders;
    /** More CSS on the host (a skin's tokens). */
    style?: string;
    /** Where the host goes (default: the body). */
    parent?: HTMLElement;
    section?: number;
    measurer?: boolean;
    /** False: letters through `fillText`, not the glyph atlas (`SurfaceOptions.atlas`). */
    atlas?: boolean;
    /** Built-in renderers replaced by name (a host's). */
    replace?: Record<string, LayerRenderer>;
  } = {},
): Promise<Opened> {
  const manual = manualScheduler(0);
  const scheduler = o.still
    ? { ...manual, requestFrame: () => 1, setTimer: () => 1, cancelFrame() {}, clearTimer() {} }
    : manual;
  const spec = DefaultSeqViewerSpec();
  const first = {
    ...(spec.registry ?? [])[0],
    renderers: (KernelRenderers as LayerRenderer[]).map((r) => o.replace?.[r.name] ?? r),
  };
  const registry = [
    o.labels === false ? { ...first, labels: undefined } : first,
    ...(spec.registry ?? []).slice(1),
    ...(o.registry ?? []),
  ];
  const ctx = createSeqViewer(
    { ...spec, registry, render: o.render },
    { scheduler, loaders: o.loaders, text: fixtureMeasurer() },
  );
  const r = await ctx.setDocument(doc);
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  const section = ctx.viewports.sections.value[o.section ?? 0];
  const host = document.createElement("div");
  host.style.cssText = `width:${o.width ?? 800}px;height:${o.height ?? 300}px;${HOST};${o.style ?? ""}`;
  (o.parent ?? document.body).appendChild(host);
  const surface = mountSection(ctx, section, host, {
    dpr: o.dpr ?? 1,
    readback: true,
    measurer: o.measurer,
    atlas: o.atlas,
  });
  await ctx.graph.update();
  const advance = (ms = 16) => manual.advance(ms);
  advance();
  surface.renderNow();
  return { ctx, surface, host, section, advance };
}

export function close(t: Opened) {
  t.surface.dispose();
  t.ctx.dispose();
  t.host.remove();
}

/** Frames run (the layout's, then the surface's), then a draw now. */
export function settle(t: Opened) {
  for (let k = 0; k < 3; k++) t.advance();
  t.surface.renderNow();
}

/** The surface's canvases, by `Pass`: main, overlay, navigator, its overlay, vertical overview, its overlay. */
export const canvasesOf = (t: Opened) => Array.from(t.surface.element.querySelectorAll("canvas"));
export const pixelsOf = (c: HTMLCanvasElement) =>
  c.width * c.height === 0
    ? new Uint8ClampedArray(0)
    : (c.getContext("2d") as CanvasRenderingContext2D).getImageData(0, 0, c.width, c.height).data.slice();

/**
 * Pixels differing between two captures, compared premultiplied (as the canvas stores them) by more than `tolerance`
 * in a channel, and the first one's position. Skia does not rasterize a glyph identically every time (its glyph
 * cache), which moves a text pixel by one unit: `tolerance` 1 for anything with text.
 */
export function diff(a: Uint8ClampedArray, b: Uint8ClampedArray, width: number, tolerance = 0) {
  if (a.length !== b.length) return { count: -1, x: -1, y: -1 };
  let n = 0;
  let at = -1;
  for (let i = 0; i < a.length; i += 4) {
    const aa = a[i + 3];
    const ba = b[i + 3];
    let d = Math.abs(aa - ba);
    for (let c = 0; c < 3; c++)
      d = Math.max(d, Math.abs(Math.round((a[i + c] * aa) / 255) - Math.round((b[i + c] * ba) / 255)));
    if (d > tolerance) {
      if (at < 0) at = i / 4;
      n++;
    }
  }
  return { count: n, x: at < 0 ? -1 : at % width, y: at < 0 ? -1 : Math.floor(at / width) };
}
