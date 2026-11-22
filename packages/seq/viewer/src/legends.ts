/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ReactiveModel, type ReadonlyState, State } from '@molstar/core/reactive';
import { type RangesLay, resolveColor, type SeqTokens, type ViewLayout, viewLegends } from '@molstar/seqx-core/render';
import { type Key, type Location, location } from '@molstar/seqx-core/query';
import type { ViewerConfig } from './config.js';
import { formatValue, type Legend, legendsEqual, legendTitle } from '@molstar/seqx-core/theme/index';

// The viewer's legends (_spec/README.md "Legends"): derived from the layout's channels when the view settles, a theme,
// scale or token table changes, or the hover moves onto another entry; never per frame, compared structurally. The
// reader's color-blind substitution (`theme.cvdSafe`) is state here: the reader's, never the spec's.

/** The entry the hover is on: a legend's key and the value (text) or number the hovered item has. */
export interface LegendHover {
  readonly key: string;
  readonly value?: string;
  readonly number?: number;
}

/**
 * Each item drawn by a legend's ranges layers (their hit marks: inside a merged run too, never items read and out of
 * view), over every strip's lays (a section of blocks has one a block) and its tracks' (the first strip's): an item
 * twice is one key and one range once the set and the Location are made.
 */
function drawnOf(layout: ViewLayout, layers: readonly string[], f: (lay: RangesLay, item: number) => void): void {
  for (const s of layout.sections)
    for (const band of [...s.tracks.map((t) => t.layers), ...s.lays])
      for (const l of band ?? []) {
        if (l.type.name !== 'ranges' || !layers.includes(l.key)) continue;
        const lay = l.lay as RangesLay;
        for (let k = 0; k < lay.hn; k++)
          if (lay.hitem[k] >= 0 && lay.painted[lay.hitem[k]] !== undefined) f(lay, lay.hitem[k]);
      }
}

/** The items an entry names in view (`LegendModel.entryItems`). */
export type EntryItems = { readonly elements: Key[]; readonly location: Location };

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

  /** The legends of a settled layout, published when they differ; and the entries whose items it draws. */
  update(layout: ViewLayout): void {
    this.derived++;
    this._layout = layout;
    this._items.clear();
    this._legends.set(viewLegends(layout).legends);
    // the entries with items drawn, read once a settled frame (never a frame, never a redraw of the part)
    const drawn = new Set<string>();
    for (const legend of this._legends.value)
      drawnOf(layout, legend.layers, (lay, j) => drawn.add(`${legend.key}\u0000${lay.painted[j]}`));
    this._drawnSet = drawn;
    this._drawn.set([...drawn].sort().join('\n'));
  }

  private _layout: ViewLayout | undefined;
  /** `entryItems` of the layout read (by legend and value; null for none): a hover again on an entry makes nothing. */
  private readonly _items = new Map<string, EntryItems | null>();
  private _drawnSet = new Set<string>();
  private readonly _drawn = new State('');
  /** The entries with items drawn at the last settled frame, as a key: a part redraws when it moves (a pan, tiles). */
  readonly drawn: ReadonlyState<string> = this._drawn.asReadonly;

  /**
   * The items in view an entry names (_spec/README.md "Legends", both ways): the drawn items of the legend's ranges
   * layers whose painted value is the entry's (what `updateHover` reads the other way), with their extents on their data's
   * systems; undefined when none is drawn. What the entry's hover lights where its address names no items it holds (a
   * source's rows, a derived table).
   */
  entryItems(layer: string, value: string): EntryItems | undefined {
    const legend = this._legends.value.find((l) => l.layers.includes(layer));
    const layout = this._layout;
    if (!legend || !layout) return undefined;
    const id = `${legend.key}\u0000${value}`;
    const known = this._items.get(id);
    if (known !== undefined) return known ?? undefined;
    const items = this.readItems(layout, legend.layers, value);
    this._items.set(id, items ?? null);
    return items;
  }

  private readItems(layout: ViewLayout, layers: readonly string[], value: string): EntryItems | undefined {
    const elements: Key[] = [];
    const parts: { system: string; ranges: number[] }[] = [];
    drawnOf(layout, layers, (lay, j) => {
      if (lay.painted[j] !== value) return;
      elements.push({ owner: lay.owner, key: lay.keys[j] });
      parts.push({ system: lay.dataSystem, ranges: [lay.starts[j], Math.max(lay.ends[j], lay.starts[j] + 1)] });
    });
    if (!elements.length) return undefined;
    return { elements, location: location({ parts }) };
  }

  /**
   * Whether an address names an entry's items (the theming behavior's, where it is listed: a layer's data by its
   * legend's field); none, only what is drawn is named (`entryItems`).
   */
  addressable: ((layer: string, value: string) => boolean) | undefined = undefined;

  /**
   * Whether an entry names items: by an address, or drawn at the last settled frame (a click selects them, a hover
   * lights them); a lookup, never a walk of the layout.
   */
  names(layer: string, value: string): boolean {
    if (this.addressable?.(layer, value) === true) return true;
    const legend = this._legends.value.find((l) => l.layers.includes(layer));
    return legend !== undefined && this._drawnSet.has(`${legend.key}\u0000${value}`);
  }

  /** The entry a hit is on: the item's painted value on its layer's legend (ranges), else none. */
  updateHover(
    layout: ViewLayout | undefined,
    hit: { readonly layer?: string; readonly element?: { readonly key: string | number } } | undefined,
  ): void {
    if (layout && layout !== this._layout) {
      this._layout = layout;
      this._items.clear();
    }
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
  /** A size key's mark: what the channel gives the value, in the layer's own units; drawn as a circle by ratio. */
  readonly size?: number;
  /** A step: its cells join in the compact form. */
  readonly step?: true;
  /** The value a click or Enter selects (swatches). */
  readonly value?: string;
  /**
   * The entry names items (`LegendModel.names`): a hover lights them, a click selects them, and the pointer says so; an
   * entry with no action (a swatch alone) is no target, the default cursor (_spec/README.md "Legends", both ways).
   */
  readonly action?: true;
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
  const names = (layer: string, value: string | undefined) =>
    value !== undefined && viewer.legends.names(layer, value) ? ({ action: true } as const) : {};
  const max = options.maxItems ?? 12;
  const out: LegendRow[] = [];
  for (const l of selectLegends(viewer.legends.legends.value, options)) {
    const e = l.entries;
    let entries: LegendEntryRow[] = [];
    let more = l.more;
    if (e.kind === 'swatches') {
      entries = e.items.map((i) => ({
        label: i.label,
        name: `${legendTitle(l)}: ${i.label}`,
        color: i.color,
        ...(i.value !== undefined ? { value: i.value } : {}),
        ...names(l.layers[0], i.value),
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
          name: `${legendTitle(l)}: ${range}${l.open ? ', in view' : ''}`,
          ramp: e.stops.map((s) => s.color),
          ends: [formatValue(lo), formatValue(hi)],
          hovered: hover?.key === l.key,
        },
      ];
    } else if (e.kind === 'steps') {
      entries = e.colors.map((c, i) => {
        const label = `${formatValue(e.edges[i])}–${formatValue(e.edges[i + 1])}`;
        return { label, name: `${legendTitle(l)}: ${label}`, color: c, step: true, hovered: false };
      });
    } else {
      entries = e.items.map((i) => ({
        label: formatValue(i.value),
        name: `${legendTitle(l)}: ${formatValue(i.value)}`,
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
      title: legendTitle(l),
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
