/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for the DOM UI's tests (UX2): a real context over a document by value (no network, no `node:fs` in the
// browser), the manual clock, and a canvas stub in place of `mountSection` (a focusable element sized like a surface,
// the section's viewport sized as the engine would). Imports no React.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { ShownIndex } from "../../src/axis.js";
import { sectionLabel } from "../../src/chrome.js";
import type { SurfaceFactory } from "../../src/dom/index.js";
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  EntryKind,
  Inputs,
  Modifier,
  manualScheduler,
  PointerType,
  SeqViewerConfig,
  type SeqViewerContext,
  type SeqViewerSpec,
} from "../../src/index.js";
import { fixtureMeasurer } from "../../src/testing/index.js";

const PROTEIN = "MSEEKLKELAEKAIRLAKEGNPEEALKLLEEALRLAKESGDPELIARVLELLAEVTKQSVHAYRLDDGKWFLQPGGRVSLEMAPTDLQAAH";
export const WIDTH = 900;
export const HEIGHT = 220;

/** Two views (chains; an alignment with rows), a description with command links, two snapshots with their own. */
export function uiDoc() {
  const b = createSeqViewSpec({
    title: "UI fixture",
    description: "The kinase and its domains, as given.\n\nA second paragraph.",
  });
  const a = b.data.sequence({ ref: "a", label: "Kinase", value: PROTEIN });
  const dom = b.data.annotation({
    ref: "dom",
    geometry: "ranges",
    target: a,
    value: [
      { start: 4, stop: 30, label: "N-lobe" },
      { start: 40, stop: 80, label: "C-lobe" },
    ],
  });
  const aln = b.data.alignment({
    ref: "aln",
    value: Array.from({ length: 6 }, (_, k) => PROTEIN.slice(k, k + 60)),
    labels: ["One", "Two", "Three", "Four", "Five", "Six"],
  });
  const chains = b.view({ ref: "chains", label: "Chains" }).block({ sequence: a });
  chains.track({ ref: "t.ruler", label: "" }).ruler();
  chains.track({ ref: "t.res", label: "Residues" }).letters();
  chains.track({ ref: "t.dom", label: "Domains" }).ranges({ data: dom });
  const msa = b.view({ ref: "msa", label: "Alignment", style: "msa" }).block({ alignment: aln });
  msa.track({ label: "" }).ruler();
  msa.track({ ref: "t.rows", label: "Rows" }).rows();
  b.snapshot({
    ref: "s1",
    title: "Domains hidden",
    description: "Only the first domain is hidden.",
    ops: [{ op: "set", ref: "t.dom", params: { hidden: true } }],
  });
  b.snapshot({ ref: "s2", title: "Unchanged", ops: [] });
  // The default Description shows plain text only (seqspec:interaction#the-default-ui): both descriptions say so.
  const doc = b.getDocument();
  return {
    ...doc,
    metadata: { ...doc.metadata, description_format: "plaintext" as const },
    snapshots: doc.snapshots?.map((s) => (s.description ? { ...s, description_format: "plaintext" as const } : s)),
  };
}

/** The canvas engine's stand-in: an element as focusable as a surface, the viewport sized. */
export const stubSurface: SurfaceFactory = (ctx, section, host) => {
  const el = document.createElement("div");
  // As a surface is (`SectionSurface`): one tab stop, a group with the role description.
  el.tabIndex = 0;
  el.setAttribute("role", "group");
  el.setAttribute("aria-roledescription", "sequence view");
  el.style.cssText = "position:relative;width:100%;height:100%;outline:none";
  host.appendChild(el);
  ctx.viewports.get(section)?.setSize(WIDTH, HEIGHT, 1);
  const named = ctx.document.shown.subscribe(() => {
    const label = sectionLabel(ctx, section);
    if (el.getAttribute("aria-label") !== label) el.setAttribute("aria-label", label);
  });
  return {
    element: el,
    dispose: () => {
      named.unsubscribe();
      el.remove();
    },
  };
};

/** The default spec: the renderers draw nothing until a canvas surface mounts (the tests mount a stub). */
/**
 * The default spec with the legend expanded (`UI.Legend`) and `full` chrome (`UI.Chrome`), as the UI tests were
 * written; their defaults are collapsed and clean.
 */
export const headlessSpec = (): SeqViewerSpec => {
  const spec = DefaultSeqViewerSpec();
  const config = spec.config ?? [];
  return {
    ...spec,
    config: [...config, [SeqViewerConfig.UI.Legend, "expanded"], [SeqViewerConfig.UI.Chrome, "full"]],
  };
};

export async function openedContext(): Promise<{ ctx: SeqViewerContext; advance(ms: number): void }> {
  const scheduler = manualScheduler();
  const ctx = createSeqViewer(headlessSpec(), { scheduler, text: fixtureMeasurer() });
  const r = await ctx.setDocument(uiDoc());
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  return { ctx, advance: (ms) => scheduler.advance(ms) };
}

/** Lets pending graph work and state deliveries finish. */
export async function flush() {
  for (let i = 0; i < 3; i++) await new Promise((r) => setTimeout(r, 0));
}

/** Pointer records into `ctx.input`, as the engine's DOM input fills them. */
export function pointer(ctx: SeqViewerContext) {
  const p = Inputs.pointer();
  return (type: keyof typeof PointerType, section: string, x: number, y: number) => {
    Object.assign(p, {
      type: PointerType[type],
      section: ctx.layout.sectionIndex(section),
      x,
      y,
      button: type === "move" ? -1 : 0,
      buttons: 0,
      modifiers: 0,
      pointerId: 0,
      time: 0,
    });
    ctx.input.pointer(p);
  };
}

/** Key records into `ctx.input`, as the engine's DOM input fills them; true: the keyboard behavior handled it. */
export function keys(ctx: SeqViewerContext) {
  const k = Inputs.key();
  return (section: string, key: string, shift = false) => {
    Object.assign(k, {
      section: ctx.layout.sectionIndex(section),
      key,
      code: "",
      modifiers: shift ? Modifier.shift : 0,
    });
    return ctx.input.key(k);
  };
}

/** Surface px of a track's entry center and of its header's middle. */
export function at(ctx: SeqViewerContext, section: string, track: string) {
  const view = ctx.graph.assembled.value?.view;
  const l = ctx.layout.get(section)?.current;
  const vp = ctx.viewports.get(section)?.live.current;
  if (!view || !l || !vp) throw new Error("not laid out");
  for (let e = 0; e < l.entryCount; e++) {
    if (l.entryKind[e] !== EntryKind.track) continue;
    if (view.set.strings.get(view.resolved.tracks.ref[l.entryTrack[e]]) !== track) continue;
    const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - vp.scrollY;
    return { y: top + l.entryHeight[e] / 2, header: l.headerWidth / 2 };
  }
  throw new Error(`no entry for ${track}`);
}
/** A track's header lines in view (a rows track's rows, else its entry): surface px of each top, height and middle. */
export function lines(ctx: SeqViewerContext, section: string, track: string) {
  const view = ctx.graph.assembled.value?.view;
  const l = ctx.layout.get(section)?.current;
  const vp = ctx.viewports.get(section)?.live.current;
  if (!view || !l || !vp) throw new Error("not laid out");
  const out: { top: number; height: number; y: number; header: number }[] = [];
  for (let e = 0; e < l.entryCount; e++) {
    if (l.entryKind[e] !== EntryKind.track) continue;
    if (view.set.strings.get(view.resolved.tracks.ref[l.entryTrack[e]]) !== track) continue;
    const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - vp.scrollY;
    const height = l.entryHeight[e];
    if (top >= 0 && top + height <= vp.height)
      out.push({ top, height, y: top + height / 2, header: l.headerWidth / 2 });
  }
  return out;
}
/** Surface px of a model index's cell center. */
export function xOf(ctx: SeqViewerContext, section: string, index: number) {
  const v = ctx.layout.axis(section);
  if (!v) throw new Error("no axis");
  const k = v.shownOf(0, index);
  return (v.x(0, k) + v.x(0, (k + 1) as ShownIndex)) / 2;
}

/**
 * Every live subscription the UIs can hold on the context: its states and events (by `subscriberCount`), the per-key
 * states in its managers' maps included.
 */
export function subscriptions(ctx: SeqViewerContext): number {
  const seen = new Set<object>();
  let n = 0;
  const walk = (o: unknown, depth: number) => {
    if (o === null || typeof o !== "object" || seen.has(o) || depth > 6 || ArrayBuffer.isView(o)) return;
    seen.add(o);
    if ("subscriberCount" in o) {
      const c = (o as { subscriberCount: unknown }).subscriberCount;
      n += typeof c === "function" ? (c as () => number).call(o) : (c as number);
      return;
    }
    if (o instanceof Map || o instanceof Set) for (const v of o.values()) walk(v, depth + 1);
    else if (Array.isArray(o)) for (const v of o.length < 256 ? o : []) walk(v, depth + 1);
    else for (const k of Object.keys(o)) walk((o as Record<string, unknown>)[k], depth + 1);
  };
  walk(ctx, 0);
  return n;
}
