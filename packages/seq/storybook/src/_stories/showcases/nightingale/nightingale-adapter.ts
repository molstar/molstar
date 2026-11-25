/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column } from '@molstar/core/repr/tree';
import type { DataStore, RangesField, Tables, Want } from '@molstar/seqx-core/data';
import type { Location } from '@molstar/seqx-core/query';
import { DarkTokens, LightTokens, layerSource } from '@molstar/seqx-core/render';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import { builtInSchemes, colorHex, compileChannel, type ThemeBody, type ThemeDatum } from '@molstar/seqx-core/theme';
import type { ResolvedLayer, ResolvedView } from '@molstar/seqx-core/view';

// Nightingale reads the resolved view (seqx-core _spec "Resolved view", A second consumer): the first section's
// letters, ranges and values come from the resolved view and the store through public entries only (`readSection`),
// then plain functions turn them into the shapes Nightingale's elements take (`toFeatures`, `toLineData`, the sequence
// string, `toHighlight`). Nothing here imports Nightingale: the shapes are declared as the elements read them, so the
// mapping runs in Node.

/** A range as the store holds it: `[start, end)`, 0-based, on the section's system, with the color its layer's theme gives. */
export interface RangeRow {
  readonly id: string;
  readonly label: string;
  readonly start: number;
  readonly end: number;
  /** `#rrggbb`. */
  readonly fill: string;
}

/** One position's value, 0-based; NaN where the data has none. */
export interface ValueRow {
  readonly position: number;
  readonly value: number;
}

/** What one section's data layers say, read from the store. */
export interface SectionRead {
  /** The section's system address (a location's part names it). */
  readonly system: string;
  readonly length: number;
  readonly letters: string;
  readonly ranges: readonly RangeRow[];
  readonly values: readonly ValueRow[];
  /** `#rrggbb`: the values layer's color. */
  readonly valueColor: string;
  /** Tiles not loaded when read. */
  readonly pending: number;
}

export type ThemeName = 'light' | 'dark';

const Schemes = builtInSchemes();
const hex = (rgba: number) => colorHex(rgba).slice(0, 7);

/** The token colors as `#rrggbb` for a theme: what Nightingale's attributes and styles are set from. */
export function themeColors(theme: ThemeName) {
  const t = (theme === 'dark' ? DarkTokens : LightTokens).colors;
  return {
    text: hex(t['msui-text']),
    muted: hex(t['msui-text-muted']),
    track: hex(t['seq-track']),
    line: hex(t['seq-value']),
    highlight: hex(t['seq-window']),
  };
}

const cell = (column: Column | undefined, i: number): string | number | null =>
  !column ? null : column instanceof CategoryColumn ? column.get(i) : column[i];

const fieldOf =
  (fields: readonly RangesField[] | undefined, row: number): ThemeDatum['field'] =>
  (name) => {
    const f = fields?.find((x) => x.name === name);
    const v = f ? cell(f.column, row) : undefined;
    return v === null ? undefined : v;
  };

async function tablesOf(store: DataStore, source: Want['source'], system: string, start: number, end: number) {
  const level = store.coverage(source)?.levels[0] ?? 1;
  const want: Want = { source, system, start, end, level };
  const wanter = store.wanter();
  wanter.want([want]);
  await wanter.whenSettled(2000);
  const tables: Tables[] = [];
  let pending = 0;
  for (const key of store.keysOf(want)) {
    const e = store.get(key);
    if (e?.status === 'loaded') tables.push(e.tables);
    else if (e?.status !== 'empty') pending++;
  }
  wanter.dispose();
  return { tables, pending };
}

function colorOf(layer: ResolvedLayer, theme: ThemeName, fallback: string) {
  return compileChannel('color', layer.theme.color as ThemeBody | undefined, fallback as unknown as ThemeBody, {
    tokens: theme === 'dark' ? DarkTokens : LightTokens,
    schemes: Schemes,
  });
}

const datum = (over: Partial<ThemeDatum>): ThemeDatum => ({
  field: () => undefined,
  value: Number.NaN,
  letter: '',
  start: 0,
  end: 0,
  ...over,
});

/**
 * Reads the first section's letters, ranges and values layers (the first of each) for the whole system, once the store
 * has them. Colors are the layers' theme forms compiled for `theme`.
 */
export async function readSection(
  model: StaticModel,
  view: ResolvedView,
  store: DataStore,
  theme: ThemeName = 'light',
): Promise<SectionRead> {
  const section = view.sections[0];
  if (!section) throw new Error('the view has no section');
  let length = 0;
  let letters = '';
  const ranges: RangeRow[] = [];
  const values: ValueRow[] = [];
  let valueColor = themeColors(theme).line;
  let pending = 0;
  for (const track of section.tracks)
    for (const layer of track.layers) {
      const access = layerSource(model, layer);
      // a layer reached through a path (a follower's) is carried by the path's maps: not read here
      if (!access || layer.path) continue;
      length = Math.max(length, access.length ?? 0);
      const read = await tablesOf(store, access.source, access.system, 0, access.length ?? 1 << 20);
      pending += read.pending;
      if (layer.type === 'letters' && !letters) {
        const chars: string[] = [];
        for (const t of read.tables) {
          const l = t.letters;
          if (!l) continue;
          for (let i = 0; i < l.letters.length; i++)
            chars[l.start + i] = typeof l.letters === 'string' ? l.letters[i] : String.fromCharCode(l.letters[i]);
        }
        letters = Array.from(chars, (c) => c ?? '-').join('');
      } else if (layer.type === 'ranges' && !ranges.length) {
        const color = colorOf(layer, theme, 'token:range');
        const seen = new Set<string | number>();
        for (const t of read.tables) {
          const r = t.ranges;
          if (r?.form !== 'list') continue;
          for (let i = 0; i < r.start.length; i++) {
            const key = r.id?.get(i) ?? r.row?.[i];
            if (key !== undefined && key !== null) {
              if (seen.has(key)) continue;
              seen.add(key);
            }
            const field = fieldOf(r.fields, i);
            const label = field('label') ?? r.id?.get(i) ?? '';
            const fill = color.constant ?? color.at(datum({ field, start: r.start[i], end: r.end[i] }));
            ranges.push({
              id: String(r.id?.get(i) ?? i),
              label: String(label),
              start: r.start[i],
              end: r.end[i],
              fill: hex(fill as number),
            });
          }
        }
        ranges.sort((a, b) => a.start - b.start || a.end - b.end);
      } else if (layer.type === 'values' && !values.length) {
        const color = colorOf(layer, theme, 'token:value');
        if (typeof color.constant === 'number') valueColor = hex(color.constant);
        for (const t of read.tables) {
          const v = t.values;
          if (v?.form !== 'positions') continue;
          for (let i = 0; i < v.value.length; i++) values.push({ position: v.start + i, value: v.value[i] });
        }
        values.sort((a, b) => a.position - b.position);
      }
    }
  return {
    system: section.axis.system,
    length: length || letters.length,
    letters,
    ranges,
    values,
    valueColor,
    pending,
  };
}

// Nightingale's shapes, as `nightingale-track` and `nightingale-linegraph-track` take them.

/** `nightingale-track`'s feature: locations are 1-based and inclusive. */
export interface NgFeature {
  readonly accession: string;
  readonly color: string;
  readonly shape: 'roundRectangle';
  readonly type: string;
  readonly tooltipContent: string;
  readonly locations: readonly { readonly fragments: readonly { readonly start: number; readonly end: number }[] }[];
}

/** `nightingale-linegraph-track`'s line. */
export interface NgLine {
  readonly name: string;
  readonly range: readonly [number, number];
  readonly color: string;
  readonly fill: string;
  readonly lineCurve: 'curveLinear';
  readonly values: readonly { readonly position: number; readonly value: number }[];
}

/** Ranges as features: `[start, end)` 0-based becomes `start + 1 … end`, each in the color its theme gave. */
export function toFeatures(ranges: readonly RangeRow[]): NgFeature[] {
  return ranges.map((r) => ({
    accession: r.id,
    color: r.fill,
    shape: 'roundRectangle',
    type: 'domain',
    tooltipContent: `${r.label}, ${r.start + 1}–${r.end}`,
    locations: [{ fragments: [{ start: r.start + 1, end: r.end }] }],
  }));
}

/** Values as one line (positions 1-based; a position with no value left out). */
export function toLineData(values: readonly ValueRow[], name: string, color: string, fill: string): NgLine[] {
  const rows = values.filter((v) => !Number.isNaN(v.value));
  const ys = rows.map((v) => v.value);
  return [
    {
      name,
      range: [Math.min(0, ...ys), Math.max(...ys)],
      color,
      fill,
      lineCurve: 'curveLinear',
      values: rows.map((v) => ({ position: v.position + 1, value: v.value })),
    },
  ];
}

/**
 * A location as Nightingale's `highlight` attribute (`start:end[,start:end]`, 1-based and inclusive) on one system:
 * the same address the viewer's hover and selection name, so no mapping of positions is made but the base.
 */
export function toHighlight(location: Location | undefined, system: string): string {
  const out: string[] = [];
  for (const part of location?.parts ?? []) {
    if (part.system !== system) continue;
    for (let i = 0; i + 1 < part.ranges.length; i += 2) out.push(`${part.ranges[i] + 1}:${part.ranges[i + 1]}`);
  }
  return out.join(',');
}
