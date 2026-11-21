/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ReactiveModel, type ReadonlyState, State } from '@molstar/core/reactive';
import { type RangesLay, resolveColor, type SeqTokens, type ViewLayout, viewLegends } from '@molstar/seqx-core/render';
import type { ViewerConfig } from './config.js';
import { formatValue, type Legend, legendsEqual } from '@molstar/seqx-core/theme/index';

// The viewer's legends (_spec/README.md "Legends"): derived from the layout's channels when the view settles, a theme,
// scale or token table changes, or the hover moves onto another entry; never per frame, compared structurally. The
// reader's color-blind substitution (`theme.cvdSafe`) is state here: the reader's, never the spec's.

/** The entry the hover is on: a legend's key and the value (text) or number the hovered item has. */
export interface LegendHover {
  readonly key: string;
  readonly value?: string;
  readonly number?: number;
}

export class LegendModel extends ReactiveModel {
  private readonly _legends = new State<readonly Legend[]>([], legendsEqual);
  private readonly _hover = new State<LegendHover | undefined>(
    undefined,
    (a, b) => a?.key === b?.key && a?.value === b?.value && a?.number === b?.number,
  );
  private readonly _cvdSafe: State<boolean>;
  /** The legends as last published: merged, in the view's order. */
  readonly legends: ReadonlyState<readonly Legend[]>;
  /** The legend entry the hovered item reads as (outlined; a ramp's value ticked). */
  readonly hover: ReadonlyState<LegendHover | undefined>;
  /** The reader's color-blind substitution: each unsafe scheme drawn as its substitute, a palette as Okabe–Ito. */
  readonly cvdSafe: ReadonlyState<boolean>;
  /** Legends' derivations done (a test reads that a frame alone does not derive them). */
  derived = 0;

  constructor(cvdSafe = false) {
    super();
    this._cvdSafe = new State(cvdSafe);
    this.legends = this._legends.asReadonly;
    this.hover = this._hover.asReadonly;
    this.cvdSafe = this._cvdSafe.asReadonly;
  }

  /** Sets the reader's substitution (a reader's setting, never written into the spec). */
  setCvdSafe(on: boolean): void {
    this._cvdSafe.set(on);
  }

  /** The legends of a settled layout, published when they differ. */
  update(layout: ViewLayout): void {
    this.derived++;
    this._legends.set(viewLegends(layout).legends);
  }

  /** The entry a hit is on: the item's painted value on its layer's legend (ranges), else none. */
  updateHover(
    layout: ViewLayout | undefined,
    hit: { readonly layer?: string; readonly element?: { readonly key: string | number } } | undefined,
  ): void {
    if (!layout || !hit?.element || hit.layer === undefined) return this._hover.set(undefined);
    const layer = hit.layer;
    const legend = this._legends.value.find((l) => l.layers.includes(layer));
    if (!legend) return this._hover.set(undefined);
    for (const s of layout.sections)
      for (const t of s.tracks)
        for (const l of t.layers) {
          if (l.key !== layer || l.type.name !== 'ranges') continue;
          const lay = l.lay as RangesLay;
          for (let j = 0; j < lay.items; j++)
            if (lay.keys[j] === hit.element.key) {
              const value = lay.painted[j];
              return this._hover.set(value !== undefined ? { key: legend.key, value } : { key: legend.key });
            }
        }
    this._hover.set(undefined);
  }
}

/** One entry of a legend as the part shows it: a label, its color (or a ramp's), the value it selects, its name. */
export interface LegendEntryRow {
  readonly label: string;
  /** The accessible name: the value and the title ("Domains: DNA-binding"); a ramp its domain and whether open. */
  readonly name: string;
  readonly color?: number;
  readonly ramp?: readonly number[];
  /** A ramp's domain ends as written (the compact form shows them beside its bar). */
  readonly ends?: readonly [string, string];
  /** A size key's mark: the radius (px) the channel gives the value, drawn as a circle. */
  readonly size?: number;
  /** A step: its cells join in the compact form. */
  readonly step?: true;
  /** The value a click or Enter selects (swatches). */
  readonly value?: string;
  readonly hovered: boolean;
}

export interface LegendRow {
  readonly key: string;
  readonly title: string;
  /** The layer whose items an entry selects (the first the legend keys). */
  readonly layer: string;
  readonly entries: readonly LegendEntryRow[];
  readonly more: number;
  /** What the legend says once: substituted colors, hashed values listed as painted. */
  readonly note?: string;
}

/** Which legends a part shows (_spec/README.md "Legends"): `of` and `omit` name addresses (a track's, a layer's, a scale's). */
export interface LegendSelection {
  /** An address or a list of them: only their legends (every legend when absent). */
  readonly of?: string | readonly string[];
  /** Addresses whose legends are left out (the host shows them elsewhere); applied after `of`. */
  readonly omit?: readonly string[];
  /** Entries a legend lists before "+n more" (12). */
  readonly maxItems?: number;
}

/** Whether a legend is keyed by an address: its scale, or a layer at or below it. */
const keyedBy = (l: Legend, a: string) => l.scale === a || l.layers.some((x) => x === a || x.startsWith(`${a}/`));

/** The legends a selection shows. */
export function selectLegends(all: readonly Legend[], { of, omit }: LegendSelection): Legend[] {
  const some = typeof of === 'string' ? [of] : of;
  return all.filter((l) => (!some || some.some((a) => keyedBy(l, a))) && !omit?.some((a) => keyedBy(l, a)));
}

/**
 * The rows the `Legend` part shows (DOM-free, the headless twin's): the legends of `of` (an address or a list) minus
 * `omit`, else all, each entry with its accessible name, the hovered one marked.
 */
export function legendRows(viewer: { readonly legends: LegendModel }, options: LegendSelection = {}): LegendRow[] {
  const hover = viewer.legends.hover.value;
  const max = options.maxItems ?? 12;
  const out: LegendRow[] = [];
  for (const l of selectLegends(viewer.legends.legends.value, options)) {
    const e = l.entries;
    let entries: LegendEntryRow[] = [];
    let more = l.more;
    if (e.kind === 'swatches') {
      entries = e.items.map((i) => ({
        label: i.label,
        name: `${l.title}: ${i.label}`,
        color: i.color,
        ...(i.value !== undefined ? { value: i.value } : {}),
        hovered: hover?.key === l.key && hover.value === (i.value ?? i.label),
      }));
      if (entries.length > max) {
        more += entries.length - max;
        entries = entries.slice(0, max);
      }
    } else if (e.kind === 'ramp') {
      const [lo, hi] = l.domain ?? [0, 1];
      const range = `${formatValue(lo)} to ${formatValue(hi)}`;
      entries = [
        {
          label: range,
          name: `${l.title}: ${range}${l.open ? ', in view' : ''}`,
          ramp: e.stops.map((s) => s.color),
          ends: [formatValue(lo), formatValue(hi)],
          hovered: hover?.key === l.key,
        },
      ];
    } else if (e.kind === 'steps') {
      entries = e.colors.map((c, i) => {
        const label = `${formatValue(e.edges[i])}–${formatValue(e.edges[i + 1])}`;
        return { label, name: `${l.title}: ${label}`, color: c, step: true, hovered: false };
      });
    } else {
      entries = e.items.map((i) => ({
        label: formatValue(i.value),
        name: `${l.title}: ${formatValue(i.value)}`,
        size: i.size,
        hovered: false,
      }));
    }
    const notes = [
      ...(l.notes.includes('substituted') ? ['colors substituted'] : []),
      ...(l.notes.includes('unlisted') ? ['values in view'] : []),
    ];
    out.push({
      key: l.key,
      title: l.title,
      layer: l.layers[0],
      entries,
      more,
      ...(notes.length ? { note: notes.join('; ') } : {}),
    });
  }
  return out;
}

/**
 * A highlight channel's color: `highlight.colors[name]`, else `highlight.color`, resolved through the tokens (the viewer
 * draws it, the Legend part shows it).
 */
export function highlightColorOf(
  viewer: { readonly config: ViewerConfig; readonly tokens: SeqTokens },
  name: string,
): number {
  const { color, colors } = viewer.config.highlight;
  return (
    resolveColor(Object.hasOwn(colors, name) ? colors[name] : color, viewer.tokens) ??
    resolveColor(color, viewer.tokens) ??
    viewer.tokens.colors['seq-cat-5']
  );
}

export type LegendChannels = 'reserve' | 'collapse';
export const LegendChannelsModes: readonly LegendChannels[] = ['reserve', 'collapse'];
/**
 * The highlight channels' row of the Legend part (DOM-free): whether it is there, for a mode and the number of active
 * channels. `reserve` has it with or without a channel (empty, it holds the line a channel's row takes: the part draws a
 * hidden copy, so the height is the real one at any text size); `collapse` only with one.
 */
export function channelsRow(mode: LegendChannels, active: number): { present: boolean } {
  return { present: mode === 'reserve' || active > 0 };
}
