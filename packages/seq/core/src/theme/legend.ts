/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashString } from '@molstar/core/data/util/hash-functions';
import { Color } from '@molstar/core/util/color/color';
import { type Legend as MolLegend, ScaleLegend, TableLegend } from '@molstar/core/util/legend';
import { toRgb } from './colors.js';
import type { Channel } from './compile.js';
import type { ChannelName } from './forms.js';
import { normalize } from './palettes.js';

// Legends (_spec/README.md "Themes", Legends): derived from a layer's resolved channel, never written; merged when they
// say the same; converted to Mol*'s `TableLegend` and `ScaleLegend` for a Mol* host's legend control.

export type LegendEntries =
  | {
      readonly kind: 'ramp';
      /** RGBA colors at 0..1 of the domain. */
      readonly stops: readonly { readonly color: number; readonly at: number }[];
      readonly ticks: readonly number[];
      /** Each tick's place on the ramp, 0..1 through the scale (and its center at 0.5 with `mid`). */
      readonly tickAt: readonly number[];
      readonly scale: 'linear' | 'log' | 'sqrt';
      readonly mid?: number;
    }
  | { readonly kind: 'steps'; readonly edges: readonly number[]; readonly colors: readonly number[] }
  | {
      readonly kind: 'swatches';
      readonly items: readonly { readonly label: string; readonly color: number; readonly value?: string }[];
    }
  | { readonly kind: 'sizes'; readonly items: readonly { readonly value: number; readonly size: number }[] };

export interface Legend {
  /** What two legends that say the same share: they merge into one. */
  readonly key: string;
  readonly title: string;
  readonly channel: ChannelName;
  /** The scale's address when the channel names a scale node. */
  readonly scale?: string;
  /** The field the channel reads. */
  readonly field?: string;
  readonly entries: LegendEntries;
  /** The domain as resolved; `open`: an end follows the data in view ("in view"). */
  readonly domain?: readonly [number, number];
  readonly open: boolean;
  /** `vocabulary`: every listed value; `painted`: what was drawn (EC-204). */
  readonly source: 'vocabulary' | 'painted';
  /** Entries past `maxItems`. */
  readonly more: number;
  /** `unlisted`: values the palette does not list take hashed colors, listed only once painted. `substituted`: a reader's color-blind substitution. */
  readonly notes: readonly ('unlisted' | 'substituted')[];
  /** The layers (keys) it keys. */
  readonly layers: readonly string[];
}

export interface LegendOptions {
  /** The legend's title: the layer's label (its track's). */
  readonly title: string;
  readonly layer: string;
  /** The values drawn (categories painted in view), in order. */
  readonly painted?: Iterable<string>;
  /** Entries listed before `more` counts the rest: the part's prop (12). */
  readonly maxItems?: number;
}

/** A log scale's ticks: its decades inside the domain (the ends when none is). */
export function logTicks(lo: number, hi: number): number[] {
  if (!(lo > 0) || !(hi > lo)) return niceTicks(lo, hi);
  const out: number[] = [];
  for (let e = Math.ceil(Math.log10(lo)); e <= Math.floor(Math.log10(hi)); e++) out.push(10 ** e);
  return out.length ? out : [lo, hi];
}

/** Nice ticks over a domain: steps of 1, 2, 5 × 10ⁿ, at most about five. */
export function niceTicks(lo: number, hi: number, count = 5): number[] {
  if (!(hi > lo)) return [lo];
  const raw = (hi - lo) / count;
  const p = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((m) => m * p).find((s) => s >= raw) ?? 10 * p;
  const out: number[] = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + step * 1e-9; t += step) out.push(Number(t.toPrecision(12)));
  return out;
}

/** A channel's legend, or undefined (a uniform channel says nothing a key would). */
export function legendOf(channel: Channel, options: LegendOptions): Legend | undefined {
  const src = channel.legend;
  if (!src) return undefined;
  const maxItems = options.maxItems ?? 12;
  const notes: ('unlisted' | 'substituted')[] = channel.scheme?.substituted ? ['substituted'] : [];
  let entries: LegendEntries;
  let source: Legend['source'] = 'vocabulary';
  let domain: [number, number] | undefined;
  let open = false;
  let more = 0;
  let scale: string | undefined;
  let field: string | undefined;
  if (src.kind === 'swatches') {
    const items = src.items.filter((i) => typeof i.value === 'number');
    more = Math.max(0, items.length - maxItems);
    entries = {
      kind: 'swatches',
      items: items.slice(0, maxItems).map((i) => ({ label: i.label, color: i.value as number })),
    };
  } else if (src.kind === 'sizes') {
    const [lo, hi] = src.domain;
    // the size channel's own scale (a radius by `sqrt`): one `normalize` for the marks and their key
    const at = (v: number) => {
      const t = normalize(v, src.type, lo, hi, undefined);
      return src.range[0] + (src.range[1] - src.range[0]) * (Number.isNaN(t) ? 0 : Math.max(0, Math.min(1, t)));
    };
    entries = { kind: 'sizes', items: [lo, (lo + hi) / 2, hi].map((v) => ({ value: v, size: at(v) })) };
    domain = [lo, hi];
    field = src.field;
    scale = src.scale;
  } else {
    const p = src.palette;
    field = src.field;
    scale = src.scale;
    open = src.open;
    if (p.kind === 'direct') return undefined;
    if (p.kind === 'categorical') {
      const listed = p.values ?? [];
      const values: string[] = [...listed];
      const seen = new Set(values);
      for (const v of options.painted ?? [])
        if (!seen.has(v)) {
          seen.add(v);
          values.push(v);
        }
      source = listed.length ? 'vocabulary' : 'painted';
      // hashed values cannot be listed before they are painted: said once where some are (EC-204)
      if (p.hashed && (!listed.length || values.length > listed.length)) notes.push('unlisted');
      more = Math.max(0, values.length - maxItems);
      entries = {
        kind: 'swatches',
        items: values.slice(0, maxItems).flatMap((v) => {
          const c = p.color(v);
          return c === undefined ? [] : [{ label: v, color: c, value: v }];
        }),
      };
    } else if (p.kind === 'binned') {
      entries = { kind: 'steps', edges: p.edges ?? [], colors: p.colors };
      domain = p.domain ? [p.domain[0], p.domain[1]] : undefined;
    } else {
      const [lo, hi] = p.domain ?? [0, 1];
      domain = [lo, hi];
      const at = p.at ?? [];
      const ticks = p.scale === 'log' ? logTicks(lo, hi) : niceTicks(lo, hi);
      if (p.mid !== undefined && !ticks.includes(p.mid)) {
        ticks.push(p.mid);
        ticks.sort((a, b) => a - b);
      }
      entries = {
        kind: 'ramp',
        stops: p.colors.map((color, i) => ({ color, at: at[i] ?? 0 })),
        ticks,
        tickAt: ticks.map((t) => normalize(t, p.scale, lo, hi, p.mid)),
        scale: p.scale,
        ...(p.mid !== undefined ? { mid: p.mid } : {}),
      };
    }
  }
  const title = (src.kind === 'palette' && src.label) || options.title;
  const key = `${channel.name}|${scale ?? ''}|${hashString(JSON.stringify([entries, domain, open])) >>> 0}`;
  return {
    key,
    title,
    channel: channel.name,
    ...(scale ? { scale } : {}),
    ...(field ? { field } : {}),
    entries,
    ...(domain ? { domain } : {}),
    open,
    source,
    more,
    notes,
    layers: [options.layer],
  };
}

/**
 * Legends with the same key merged into one (EC-132): its layers all of theirs, titled by what their titles share (the
 * words they start with alike), else by the first; in the order first met.
 */
export function mergeLegends(legends: readonly Legend[]): Legend[] {
  const byKey = new Map<string, Legend[]>();
  for (const l of legends) {
    const g = byKey.get(l.key);
    if (g) g.push(l);
    else byKey.set(l.key, [l]);
  }
  const out: Legend[] = [];
  for (const g of byKey.values()) {
    if (g.length === 1) {
      out.push(g[0]);
      continue;
    }
    const titles = g.map((l) => l.title.split(/\s+/));
    const shared: string[] = [];
    for (let i = 0; i < titles[0].length && titles.every((t) => t[i] === titles[0][i]); i++) shared.push(titles[0][i]);
    out.push({ ...g[0], title: shared.join(' ') || g[0].title, layers: g.flatMap((l) => l.layers) });
  }
  return out;
}

/** A number as a legend writes it: up to 4 significant digits. */
export function formatValue(v: number): string {
  return Number.isInteger(v) ? String(v) : String(Number(v.toPrecision(4)));
}

/** The legend as Mol*'s: a ramp as a `ScaleLegend`, steps and swatches as a `TableLegend`; sizes have none. */
export function toMolLegend(legend: Legend): MolLegend | undefined {
  const e = legend.entries;
  if (e.kind === 'ramp') {
    const [lo, hi] = legend.domain ?? [0, 1];
    return ScaleLegend(
      formatValue(lo),
      formatValue(hi),
      e.stops.map((s) => [Color(toRgb(s.color)), s.at] as [Color, number]),
    );
  }
  if (e.kind === 'steps')
    return TableLegend(
      e.colors.map(
        (c, i) => [`${formatValue(e.edges[i])}–${formatValue(e.edges[i + 1])}`, toRgb(c)] as [string, Color],
      ),
    );
  if (e.kind === 'swatches') return TableLegend(e.items.map((i) => [i.label, toRgb(i.color)] as [string, Color]));
  return undefined;
}

/** Two legend lists equal in what they show: the viewer publishes legends on change only. */
export function legendsEqual(a: readonly Legend[], b: readonly Legend[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++)
    if (
      a[i].key !== b[i].key ||
      a[i].title !== b[i].title ||
      a[i].more !== b[i].more ||
      a[i].layers.join() !== b[i].layers.join()
    )
      return false;
  return true;
}
