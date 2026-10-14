/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The legend's content per track (seqspec:interaction#the-legend): what each shown track's layers' themes list
// (renderer `legend`), read once per assembled view and shared by the legend's groups (`legendGroups`, chrome.ts) and a
// track header's tooltip (`headerLegend`, through the hover's label). Edge-rate: computed when the view changes, never
// per frame or hover.

import { NO_COLOR } from "./draw/backend.js";
import type { ViewModel } from "./graph.js";
import type { SeqViewerRegistry } from "./registry.js";
import type { LegendSink } from "./themes.js";

export interface LegendItem {
  readonly label: string;
  /** `#rrggbb`; null: the ink. */
  readonly color: string | null;
  /** A `GlyphId`; `NO_GLYPH`: a swatch. */
  readonly glyph: number;
}
/** A continuous palette over [lo, hi] (`LegendSink.ramp`): `stops` at `t` in 0..1, colors `#rrggbb`. */
export interface LegendRamp {
  readonly stops: readonly { readonly t: number; readonly color: string }[];
  readonly lo: number;
  readonly hi: number;
  /** What the ramp measures ("contacts (number of kinds)"); absent: nothing is said. */
  readonly caption?: string;
}
/**
 * A sample size (`LegendSink.size`): `unit` a `SizeUnit`; `channel` how it is drawn (`width` a stroke, `radius` a dot,
 * a band's share a bar); `field` the field a scale's samples are of (null: a rule's); `color` the color a ramp over the
 * same field in the same layer gives its value (null: the ink).
 */
export interface LegendSize {
  readonly label: string;
  readonly value: number;
  readonly unit: number;
  readonly channel: string | null;
  readonly field: string | null;
  readonly color: string | null;
}
/** One track's legend: its shown layers' items (each once), ramps and sizes; `layers`, those that listed any. */
export interface TrackLegend {
  readonly items: readonly LegendItem[];
  readonly ramps: readonly LegendRamp[];
  readonly sizes: readonly LegendSize[];
  readonly layers: readonly number[];
}
/** A track header's tooltip legend: at most `HEADER_LEGEND_ITEMS` items, the rest a count. */
export interface HeaderLegend {
  readonly items: readonly LegendItem[];
  readonly ramps: readonly LegendRamp[];
  readonly more: number;
}
export const HEADER_LEGEND_ITEMS = 6;

const hex = (c: number) => `#${c.toString(16).padStart(6, "0")}`;
const legends = new WeakMap<ViewModel, ReadonlyMap<number, TrackLegend>>();
const headers = new WeakMap<ViewModel, Map<number, HeaderLegend | null>>();

/** A ramp's color at t (clamped to 0..1): the nearest stop's. */
function rampColor(stops: Float64Array, t: number): string | null {
  const n = stops.length >> 1;
  if (n === 0 || Number.isNaN(t)) return null;
  const c = stops[2 * Math.round(Math.min(1, Math.max(0, t)) * (n - 1)) + 1];
  return c === NO_COLOR ? null : hex(c);
}

/** Each shown track's legend, by its resolved index, in track order; a track with nothing to list is absent. */
export function trackLegends(registry: SeqViewerRegistry, view: ViewModel): ReadonlyMap<number, TrackLegend> {
  const known = legends.get(view);
  if (known) return known;
  const { resolved: rv, prepared } = view;
  const all = new Map<number, { items: LegendItem[]; ramps: LegendRamp[]; sizes: LegendSize[]; layers: number[] }>();
  for (let i = 0; i < rv.layers.count; i++) {
    const t = rv.layers.track[i];
    const r = registry.get("renderers", prepared[i]?.renderer ?? "");
    if (rv.tracks.hidden[t] || rv.sections.hidden[rv.tracks.section[t]] || !r?.legend) continue;
    const g = all.get(t) ?? { items: [], ramps: [], sizes: [], layers: [] };
    all.set(t, g);
    const before = g.items.length + g.ramps.length + g.sizes.length;
    // This layer's ramps by field: a size sample of the same field takes the ramp's color at its value (merged).
    const byField = new Map<string, { stops: Float64Array; lo: number; hi: number }>();
    const sink: LegendSink = {
      item: (label, color, glyph) => {
        const c = color === NO_COLOR ? null : hex(color);
        if (!g.items.some((x) => x.label === label && x.color === c && x.glyph === glyph))
          g.items.push({ label, color: c, glyph });
      },
      ramp: (stops, lo, hi, caption, field) => {
        if (field !== undefined) byField.set(field, { stops, lo, hi });
        g.ramps.push({
          stops: Array.from({ length: stops.length >> 1 }, (_, k) => ({
            t: stops[2 * k],
            color: hex(stops[2 * k + 1]),
          })),
          lo,
          hi,
          ...(caption === undefined ? {} : { caption }),
        });
      },
      size: (label, value, unit, of) => {
        const r = of?.field === undefined ? undefined : byField.get(of.field);
        const color =
          r && of?.at !== undefined ? rampColor(r.stops, r.hi > r.lo ? (of.at - r.lo) / (r.hi - r.lo) : 0) : null;
        const x = { label, value, unit, channel: of?.channel ?? null, field: of?.field ?? null, color };
        // Listed once (two layers of a track over one scale: a sashimi's strands).
        if (
          !g.sizes.some((y) => y.label === label && y.value === value && y.channel === x.channel && y.field === x.field)
        )
          g.sizes.push(x);
      },
    };
    r.legend(rv.layers.params[i], prepared[i].themes, sink, prepared[i].state);
    if (g.items.length + g.ramps.length + g.sizes.length > before) g.layers.push(i);
  }
  for (const [t, g] of all) if (g.layers.length === 0) all.delete(t);
  legends.set(view, all);
  return all;
}

/** A track's legend for its header's tooltip (null: nothing to list), made once per track and view. */
export function headerLegend(registry: SeqViewerRegistry, view: ViewModel, track: number): HeaderLegend | null {
  const byTrack = headers.get(view) ?? new Map<number, HeaderLegend | null>();
  headers.set(view, byTrack);
  const known = byTrack.get(track);
  if (known !== undefined) return known;
  const g = trackLegends(registry, view).get(track);
  const out =
    g && g.items.length + g.ramps.length > 0
      ? {
          items: g.items.slice(0, HEADER_LEGEND_ITEMS),
          ramps: g.ramps,
          more: Math.max(0, g.items.length - HEADER_LEGEND_ITEMS),
        }
      : null;
  byTrack.set(track, out);
  return out;
}
