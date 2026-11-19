/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ReactiveModel, type ReadonlyState, State } from '@molstar/core/reactive';
import { type RangesLay, type ViewLayout, viewLegends } from '@molstar/seqx-core/render';
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

/**
 * The rows the `Legend` part shows (DOM-free, the headless twin's): the legends of `of` (a track's or a scale's address)
 * or all, each entry with its accessible name, the hovered one marked.
 */
export function legendRows(
  viewer: { readonly legends: LegendModel },
  options: { readonly of?: string; readonly maxItems?: number } = {},
): LegendRow[] {
  const hover = viewer.legends.hover.value;
  const max = options.maxItems ?? 12;
  const out: LegendRow[] = [];
  for (const l of viewer.legends.legends.value) {
    const of = options.of;
    if (of !== undefined && l.scale !== of && !l.layers.some((a) => a === of || a.startsWith(`${of}/`))) continue;
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
          hovered: hover?.key === l.key,
        },
      ];
    } else if (e.kind === 'steps') {
      entries = e.colors.map((c, i) => {
        const label = `${formatValue(e.edges[i])}–${formatValue(e.edges[i + 1])}`;
        return { label, name: `${l.title}: ${label}`, color: c, hovered: false };
      });
    } else {
      entries = e.items.map((i) => ({
        label: formatValue(i.value),
        name: `${l.title}: ${formatValue(i.value)}`,
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
