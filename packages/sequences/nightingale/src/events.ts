/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What `mount` decides from Nightingale's events, platform-free so tests reach it without a DOM or Nightingale: the
// axis position and ranges an event names, the arbiter between a component's item hover and the pointer's by-position
// hover, the pointer's hover deduped before it is located, the per-manager highlight dedupe (`lit`), and the merge of
// fetched values into a series.

import { type AxisMeta, AxisRuns, type Located, locate } from "./axis.ts";
import { type NightingaleAxis, positionAt } from "./labels.ts";

export type Detail = Record<string, unknown> & { feature?: Record<string, unknown> | null };

/** The first value of an object's own enumerable keys (the linegraph's feature has its series' one key). */
function onlyValue(o: Record<string, unknown>): unknown {
  for (const k in o) return o[k];
  return undefined;
}
/**
 * The 1-based axis position an event names: its highlight, a residue's position, or a line's point (`series`: the
 * linegraph's name, the one key of its feature).
 */
export function positionOf(detail: Detail, series?: string): number {
  const h = (detail.highlight ?? (detail.type === "highlight" ? detail.value : null)) as string | null;
  if (typeof h === "string" && h) return Number.parseInt(h, 10);
  const f = detail.feature;
  if (typeof f?.position === "number") return f.position;
  if (typeof f?.start === "number") return f.start;
  const line = (f && (series !== undefined ? f[series] : onlyValue(f))) as { position?: number } | undefined;
  return typeof line?.position === "number" ? line.position : -1;
}

interface Fragments {
  locations?: { fragments: { start: number; end: number }[] }[];
}
/**
 * The axis runs an event covers, into `out`: an item's fragments (Nightingale's own highlight is their hull) read from
 * our feature's numbers, else its highlight (a contact's ends, a residue) parsed by a char-code scan, else the one
 * position.
 */
export function rangesOf(detail: Detail, at: number, out: AxisRuns): AxisRuns {
  out.count = 0;
  out.sorted = true;
  const f = detail.feature;
  const locations = f?.svs ? (f as Fragments).locations : null;
  if (locations?.length) {
    for (const l of locations) for (const x of l.fragments) AxisRuns.add(out, x.start - 1, x.end);
    return AxisRuns.finish(out);
  }
  const h = detail.highlight ?? (detail.type === "highlight" ? detail.value : null);
  if (typeof h === "string" && h) return AxisRuns.parse(h, out);
  if (at > 0) AxisRuns.add(out, at - 1, at);
  return out;
}

/** Which of a component's item and the pointer owns hover: set by an item's `hover` and `leave`, never by a click. */
export interface Arbiter {
  onItem: boolean;
}
/**
 * Whether a component's event goes on to the host: a hover only for an item (a feature, a variant, a contact: the pointer
 * hovers by position otherwise), a leave only after an item's hover. A click changes neither.
 */
export function arbitrate(a: Arbiter, kind: "hover" | "click" | "leave", item: boolean): boolean {
  if (kind === "hover" && !item) return false;
  if (kind === "leave" && !a.onItem) return false;
  if (kind !== "click") a.onItem = kind === "hover";
  return true;
}

/** The pointer's hover: the manager, axis position and layer element last reported (compared before anything is located). */
export interface HoverTarget {
  manager: unknown;
  position: number;
  layer: unknown;
}
export const HoverResult = { off: -1, same: 0, moved: 1 } as const;
/**
 * The pointer at `x` in a row's data cell, whose elements show `view`: `moved` (a hover change, `out` written),
 * `same` (the position, manager and layer element as last reported: nothing is located) or `off` (a gutter, or no
 * shown position there). Allocates nothing.
 */
export function hoverAt(
  state: HoverTarget,
  manager: unknown,
  meta: AxisMeta,
  view: NightingaleAxis,
  x: number,
  layer: unknown,
  out: Located,
): number {
  const at = positionAt(view, x);
  if (at < 0) return HoverResult.off;
  if (state.manager === manager && state.position === at && state.layer === layer) return HoverResult.same;
  if (!locate(meta, at - 1, out)) return HoverResult.off;
  state.manager = manager;
  state.position = at;
  state.layer = layer;
  return HoverResult.moved;
}

/**
 * What one manager shows: the runs last sent (`last`) and its `highlight` as Nightingale reports it (`lit`, its own hover
 * included, `own`: it came from Nightingale, so `last` is not what is shown).
 */
export interface Light {
  readonly last: AxisRuns;
  lit: string | null;
  own: boolean;
}
export const newLight = (): Light => ({ last: AxisRuns.create(), lit: null, own: false });
/** Nightingale reported this `highlight` for the manager (a child's hover, a residue). */
export function noteLit(l: Light, highlight: string | null | undefined): void {
  l.lit = highlight || null;
  l.own = true;
}
/**
 * The highlight to send a manager for `runs`: `undefined` when nothing changes (the runs equal the last sent and are what
 * is shown, or Nightingale already shows the string), else the string (null: clear). Builds a string only for a changed
 * manager: each change re-applies the manager's attributes to every element.
 */
export function nextHighlight(l: Light, runs: AxisRuns): string | null | undefined {
  if (!l.own && AxisRuns.equal(l.last, runs)) return undefined;
  const text = AxisRuns.text(runs) || null;
  AxisRuns.copy(l.last, runs);
  const same = text === l.lit;
  l.lit = text;
  l.own = false;
  return same ? undefined : text;
}

/** A linegraph's series, filled only where it has no value: fetched values never replace inline ones. */
export function mergeSeries(
  series: readonly { position: number; value: number | null }[],
  fetched: readonly (number | null)[],
): { position: number; value: number | null }[] {
  return series.map((p, k) => (p.value !== null ? p : { position: p.position, value: fetched[k] ?? null }));
}
