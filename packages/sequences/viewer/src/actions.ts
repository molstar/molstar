/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Actions (seqspec:viewer#actions): what the UI draws as menu items and switches and what hosts map into their own
// command systems. Each is a patch proposed to the document port, or a runtime change; `ctx.actions` runs them by name,
// with the origin that triggered them. The document port routes each patch. The runtime ones change `ctx.ui`, the
// viewports or the streams (`section.style`, `group.toggle`, `clade.toggle`, `viewport.fit`, `stream.retry`).

import { canonical, type DocNode, type Op, type SeqViewSpecDoc, type ThemeRef } from "@molstar/seq-core/format";
import { layerSchemes, overrideTheme } from "./color-schemes.js";
import { SeqViewerConfig } from "./config.js";
import type { SeqViewerContext } from "./context.js";
import type { Documents } from "./document.js";
import type { Origin } from "./state.js";

export interface ActionEntry<T = unknown> {
  readonly name: string;
  readonly title: string;
  readonly key?: string;
  when(ctx: SeqViewerContext, target: T): boolean;
  /** Its title for this target where it depends on the state (the hover card's button); absent: `title`. */
  titleFor?(ctx: SeqViewerContext, target: T): string;
  /** Usually `ctx.document.propose(patch, title, origin)`. */
  run(ctx: SeqViewerContext, target: T, origin: Origin): void;
}
export interface ActionManager {
  /** False: no such action, or its `when` is false. `origin`: "api" by default. */
  run(name: string, target: unknown, origin?: Origin): boolean;
  can(name: string, target: unknown): boolean;
}

export function createActions(ctx: SeqViewerContext): ActionManager {
  const entry = (name: string, target: unknown) => {
    const a = ctx.registry.get("actions", name) as ActionEntry | undefined;
    return a?.when(ctx, target) ? a : null;
  };
  return {
    can: (name, target) => entry(name, target) !== null,
    run(name, target, origin = "api") {
      const a = entry(name, target);
      a?.run(ctx, target, origin);
      return a !== null;
    },
  };
}

/** The node and its ancestors, root first; null: no such ref. */
export function nodePath(doc: SeqViewSpecDoc | null, ref: string): DocNode[] | null {
  const walk = (n: DocNode, path: DocNode[]): DocNode[] | null => {
    if (n.ref === ref) return [...path, n];
    for (const c of n.children ?? []) {
      const r = walk(c, [...path, n]);
      if (r) return r;
    }
    return null;
  };
  return doc ? walk(doc.root, []) : null;
}
const nodeOf = (ctx: SeqViewerContext, ref: unknown, kind: string) =>
  typeof ref === "string" ? (nodePath(ctx.document.shown.value, ref)?.at(-1)?.kind === kind ? ref : null) : null;
const propose = (ctx: SeqViewerContext, title: string, ops: Op[], origin: Origin) =>
  void ctx.document.propose({ format: "seq_view_spec_patch", ops }, title, origin);

/** A track's layers' data refs. */
export function trackDataRefs(track: DocNode): string[] {
  const out: string[] = [];
  for (const l of track.children ?? [])
    if (l.kind === "layer" && l.params?.type !== "rows") out.push(...[l.params?.data ?? []].flat().map(String));
  return out;
}
/**
 * What a track's layers show, as `positions` entries: each data ref as it is, or with the layer's `filter` as
 * `{ data, filter }` (a track of "Motif" over a shared feature set covers those rows, not the set), each once (a
 * filter compared as written, its keys in any order).
 */
export function trackPositions(track: DocNode): (string | { data: string; filter: unknown })[] {
  const out: (string | { data: string; filter: unknown })[] = [];
  const seen = new Set<string>();
  for (const l of track.children ?? []) {
    if (l.kind !== "layer" || l.params?.type === "rows") continue;
    const filter = l.params?.filter;
    for (const data of [l.params?.data ?? []].flat().map(String)) {
      const key = canonical([data, filter ?? null]);
      if (seen.has(key)) continue;
      seen.add(key);
      out.push(filter ? { data, filter } : data);
    }
  }
  return out;
}
/** A track's layer refs. */
const layersOf = (ctx: SeqViewerContext, track: string): string[] =>
  typeof track !== "string"
    ? []
    : (nodePath(ctx.document.shown.value, track)
        ?.at(-1)
        ?.children?.flatMap((l) => (l.kind === "layer" && l.ref ? [l.ref] : [])) ?? []);
const sectionOf = (ctx: SeqViewerContext, track: string) =>
  nodePath(ctx.document.shown.value, track)?.findLast((n) => n.kind === "section")?.ref ?? null;

type SectionValue = { readonly section: string; readonly value: unknown };
const isSectionValue = (ctx: SeqViewerContext, t: unknown): t is SectionValue =>
  t !== null && typeof t === "object" && nodeOf(ctx, (t as SectionValue).section, "section") !== null;

function action<T>(name: string, title: string, when: ActionEntry<T>["when"], run: ActionEntry<T>["run"]) {
  return { name, title, when, run } as ActionEntry;
}
const trackSet = (name: string, title: string, hidden: boolean) =>
  action<string>(
    name,
    title,
    (ctx, t) => nodeOf(ctx, t, "track") !== null,
    (ctx, t, o) => propose(ctx, title, [{ op: "set", ref: t, params: { hidden } }], o),
  );
const sectionParam = (name: string, title: string, param: string) =>
  action<SectionValue>(name, title, isSectionValue, (ctx, t, o) =>
    propose(ctx, title, [{ op: "set", ref: t.section, params: { [param]: t.value } }], o),
  );
/** `layer.colorScheme`'s target: a scheme that fits the layer, or null (its own color); `write`: into the document. */
export type ColorSchemeTarget = { readonly layer: string; readonly scheme: string | null; readonly write?: boolean };
const isSchemeTarget = (ctx: SeqViewerContext, t: unknown): t is ColorSchemeTarget => {
  const x = t as ColorSchemeTarget | null;
  if (x === null || typeof x !== "object" || nodeOf(ctx, x.layer, "layer") === null) return false;
  return x.scheme === null || layerSchemes(ctx, x.layer).some((e) => e.name === x.scheme);
};
/** `layer.colorRange`'s target: the layer's color range at runtime, `lo < hi` (any, not the data's); null: its own. */
export type ColorRangeTarget = { readonly layer: string; readonly range: readonly [number, number] | null };
const isRangeTarget = (ctx: SeqViewerContext, t: unknown): t is ColorRangeTarget => {
  const x = t as ColorRangeTarget | null;
  if (x === null || typeof x !== "object" || nodeOf(ctx, x.layer, "layer") === null) return false;
  return x.range === null || (x.range.length === 2 && Number.isFinite(x.range[0]) && x.range[0] < x.range[1]);
};
/** `view.legend`'s target: the legend open, closed, or the other way. */
export type LegendTarget = { readonly open: boolean | "toggle" };
const isLegendTarget = (_ctx: SeqViewerContext, t: unknown): t is LegendTarget => {
  const o = (t as LegendTarget | null)?.open;
  return typeof o === "boolean" || o === "toggle";
};
const always = () => true;
const history = (ctx: SeqViewerContext) => (ctx.document as Documents).port.history;

export const DefaultActions: readonly ActionEntry[] = [
  trackSet("track.hide", "Hide track", true),
  trackSet("track.show", "Show track", false),
  action<string>(
    "track.remove",
    "Remove track",
    (ctx, t) => nodeOf(ctx, t, "track") !== null,
    (ctx, t, o) => propose(ctx, "Remove track", [{ op: "remove", ref: t }], o),
  ),
  action<string>(
    "track.filterView",
    "Show only what this track covers",
    (ctx, t) => {
      const track = typeof t === "string" ? nodePath(ctx.document.shown.value, t)?.at(-1) : null;
      return track?.kind === "track" && trackDataRefs(track).length > 0 && sectionOf(ctx, t) !== null;
    },
    (ctx, t, o) => {
      const refs = trackPositions(nodePath(ctx.document.shown.value, t)?.at(-1) as DocNode);
      const positions = refs.length === 1 ? refs[0] : refs;
      propose(ctx, "Filter view", [{ op: "set", ref: sectionOf(ctx, t) as string, params: { positions } }], o);
    },
  ),
  action<string>(
    "section.clearFilter",
    "Show all positions",
    (ctx, t) => nodeOf(ctx, t, "section") !== null,
    (ctx, t, o) => propose(ctx, "Clear filter", [{ op: "set", ref: t, params: { positions: null } }], o),
  ),
  sectionParam("section.numbering", "Numbering", "numbering"),
  sectionParam("section.copies", "Copies", "copies"),
  sectionParam("section.ligandForm", "Ligand form", "ligand_form"),
  sectionParam("section.rowSizing", "Row sizing", "row_sizing"),
  action<SectionValue>("section.style", "Style", isSectionValue, (ctx, t) =>
    ctx.ui.setStyle(t.section, t.value === null ? null : String(t.value)),
  ),
  action<ColorSchemeTarget>("layer.colorScheme", "Color scheme", isSchemeTarget, (ctx, t, o) => {
    if (!t.write || t.scheme === null) return ctx.ui.setColorScheme(t.layer, t.scheme);
    const own = (nodePath(ctx.document.shown.value, t.layer)?.at(-1)?.params?.color ?? null) as ThemeRef | null;
    ctx.ui.setColorScheme(t.layer, null);
    propose(ctx, "Color scheme", [{ op: "set", ref: t.layer, params: { color: overrideTheme(own, t.scheme) } }], o);
  }),
  action<ColorRangeTarget>("layer.colorRange", "Color range", isRangeTarget, (ctx, t) =>
    ctx.ui.setColorRange(t.layer, t.range),
  ),
  // The legend open or closed (seqspec:interaction#the-legend): UI state, kept across patches; "toggle" from what is
  // shown now.
  action<LegendTarget>("view.legend", "Legend", isLegendTarget, (ctx, t) =>
    ctx.ui.setLegendExpanded(
      t.open === "toggle"
        ? !(ctx.ui.legendExpanded.value ?? ctx.config.get(SeqViewerConfig.UI.Legend) === "expanded")
        : t.open,
    ),
  ),
  action<string>(
    "view.switch",
    "Switch view",
    (ctx, t) => ctx.views.list.value.some((v) => v.ref === t),
    (ctx, t, o) => ctx.views.setCurrent(t, o),
  ),
  action<string>(
    "group.toggle",
    "Expand or collapse",
    // A group, or a track with a split layer (its rows, seqspec:format#layer-types) or a collapsing one (its members).
    (ctx, t) =>
      nodeOf(ctx, t, "group") !== null ||
      (nodeOf(ctx, t, "track") !== null &&
        !!nodePath(ctx.document.shown.value, t)
          ?.at(-1)
          ?.children?.some((l) => l.params?.split != null || l.params?.collapse != null)),
    (ctx, t) => ctx.ui.toggleGroup(t),
  ),
  action<{ alignment: string; node: number }>(
    "clade.toggle",
    "Expand or collapse clade",
    (_ctx, t) => typeof t?.alignment === "string" && t.node >= 0,
    (ctx, t) => ctx.ui.toggleClade(t.alignment, t.node),
  ),
  action("selection.clear", "Clear selection", always, (ctx, _t, o) => ctx.selection.clear(o)),
  action<string | undefined>("viewport.fit", "Fit", always, (ctx, t) => {
    // biome-ignore lint/suspicious/noFocusedTests: the viewport's verb, not a focused test.
    for (const ref of t ? [t] : ctx.viewports.sections.value) ctx.viewports.get(ref)?.fit({});
  }),
  action<string>(
    "stream.retry",
    "Retry",
    (ctx, t) => layersOf(ctx, t).some((l) => ctx.streams.hasFailed(l)),
    (ctx, t) => {
      for (const l of layersOf(ctx, t)) ctx.streams.retry(l);
    },
  ),
  action(
    "history.undo",
    "Undo",
    (ctx) => history(ctx) && ctx.history.canUndo.value,
    (ctx) => void ctx.history.undo(),
  ),
  action(
    "history.redo",
    "Redo",
    (ctx) => history(ctx) && ctx.history.canRedo.value,
    (ctx) => void ctx.history.redo(),
  ),
];
