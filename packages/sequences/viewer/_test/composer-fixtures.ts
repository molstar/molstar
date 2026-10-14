/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for the composer (CS1): a context showing a document with the kernel renderers (draw-fixtures' `drawn`),
// its composer, recording backends with the capabilities a test needs (the recording backend has none: whole frames),
// and readers of a recording: ops by code, texts, rects in a box.

import type { Loaders, SeqViewSpecDoc } from "@molstar/seq-core/format";
import { BackendCaps, Dirty, type RenderProps, SectionComposer } from "../src/draw/index.ts";
import { KernelRenderers, type LayerRenderer, type SeqViewerRegistryEntry } from "../src/index.ts";
import { Op, RecordingBackend } from "../src/testing/index.ts";
import { opened } from "./layout-fixtures.ts";

/** A recording backend that reports `caps` (copies and back targets logged, not applied). */
export function recording(caps = 0): RecordingBackend {
  const r = new RecordingBackend();
  Object.defineProperty(r, "caps", { value: caps });
  return r;
}
export const COPY_BACK = BackendCaps.Copy | BackendCaps.Back;

export interface ComposedOptions {
  width?: number;
  height?: number;
  dpr?: number;
  render?: Partial<RenderProps>;
  section?: string;
  /** Renderers by name in place of the kernel's (`null`: left out), and more registry entries. */
  replace?: Record<string, LayerRenderer | null>;
  registry?: SeqViewerRegistryEntry[];
  loaders?: Loaders;
}

/** A context showing `doc` (the kernel renderers), a section's composer (default the first), and its frames. */
export async function composed(doc: SeqViewSpecDoc, o: ComposedOptions = {}) {
  const renderers = (KernelRenderers as LayerRenderer[])
    .map((r) => (o.replace && r.name in o.replace ? o.replace[r.name] : r))
    .filter((r): r is LayerRenderer => r !== null);
  const t = await opened(doc, { renderers, render: o.render, registry: o.registry, loaders: o.loaders });
  for (const ref of t.ctx.viewports.sections.value)
    t.ctx.viewports.get(ref)?.setSize(o.width ?? 1000, o.height ?? 400, o.dpr ?? 1);
  await t.ctx.graph.update();
  const section = o.section ?? t.ctx.viewports.sections.value[0];
  const composer = new SectionComposer(t.ctx, section);
  /** A frame of `flags` into `into` (reset first). */
  const frame = (into = recording(), flags: number = Dirty.tracks) => {
    into.reset();
    composer.frame(into, flags);
    return into;
  };
  return { ...t, section, composer, frame, vp: t.ctx.viewports.get(section) };
}

/** The ops of a recording as their codes' names, in order. */
export function opNames(r: RecordingBackend): string[] {
  const names = Object.fromEntries(Object.entries(Op).map(([k, v]) => [v, k]));
  return Array.from({ length: r.opCount }, (_, k) => names[r.ops[4 * k]]);
}
/** An op's arguments (k: its index in the log). */
export function opArgs(r: RecordingBackend, k: number, n: number): number[] {
  const a = r.ops[4 * k + 3];
  return Array.from(r.args.subarray(a, a + n));
}
/** The texts of a recording, in order. */
export function textsOf(r: RecordingBackend): string[] {
  return r.textItems.text.slice(0, r.textItems.n);
}
