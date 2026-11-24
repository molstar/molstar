/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RangesField, Tables } from '@molstar/seqx-core/data';
import type { DataStore } from '@molstar/seqx-core/data';
import { LightTokens, layerSource } from '@molstar/seqx-core/render';
import { complementCode, type StaticModel } from '@molstar/seqx-core/resolve';
import {
  builtInSchemes,
  type Channel,
  colorHex,
  compileChannel,
  lettersDefault,
  type ThemeBody,
  type ThemeDatum,
} from '@molstar/seqx-core/theme';
import {
  cellOffset,
  nearestOffsetOf,
  type ResolvedLayer,
  type ResolvedSection,
  type ResolvedView,
} from '@molstar/seqx-core/view';
import { cell, readTables } from './read.js';

// The second consumer of the resolved view (EC-183; `_spec` Resolved view): letters, ranges, values and variants of one
// section, as plain SVG, from the resolved view, the static model and the store, through the public entries of
// `@molstar/seqx-core` and nothing else. Its own tiny layout (a label column, a plot, tracks stacked) and no draw list
// of core's: what it draws is a list of marks, each with the positions it stands for, the x it lands at, its label and
// its fill, so a test reads what is where.

const Schemes = builtInSchemes();

const LabelWidth = 96;
const Gap = 4;
const LaneHeight = 14;
/** A cell this wide (px) shows its letter. */
const LetterCell = 12;

export type MarkKind = 'letter' | 'range' | 'value' | 'variant';

/** One thing drawn: the positions it stands for (`[position, end)` of its data's system) and where it landed. */
export interface Mark {
  readonly kind: MarkKind;
  /** Its layer's address. */
  readonly layer: string;
  readonly position: number;
  readonly end: number;
  readonly x: number;
  readonly width: number;
  /** `#rrggbb`. */
  readonly fill: string;
  /** A letter, a range's label (else its id). */
  readonly text?: string;
  readonly value?: number;
  /** A bar's height as a share of its track's. */
  readonly share?: number;
}

export interface Drawn {
  readonly svg: string;
  readonly marks: readonly Mark[];
  /** The plot's left edge and width, px: the window `[start, end)` (positions) spans it. */
  readonly plot: { readonly x: number; readonly width: number };
  /** Tiles not loaded when drawn: such a layer is drawn as far as it is read. */
  readonly pending: number;
  /** Layers this consumer does not draw (its subset: letters, ranges, values, variants), by address. */
  readonly skipped: readonly string[];
}

export interface DrawOptions {
  readonly width: number;
  /** The positions shown, `[start, end)` of the section's system. */
  readonly window: readonly [number, number];
  /** The section's address (the view's first when absent). */
  readonly section?: string;
}

const Kinds = new Set(['letters', 'ranges', 'values', 'variants']);

const hex = (rgba: number | string | undefined, fallback = '#000000') =>
  typeof rgba === 'number' ? colorHex(rgba).slice(0, 7) : fallback;

/** A layer's channel: its theme form as the resolved view carries it, compiled by the theme entry with the type's default. */
function channelOf(
  layer: ResolvedLayer,
  name: 'color',
  fallback: ThemeBody,
  alphabet?: 'protein' | 'dna' | 'rna',
): Channel {
  return compileChannel(name, layer.theme[name] as ThemeBody | undefined, fallback, {
    tokens: LightTokens,
    schemes: Schemes,
    ...(alphabet ? { alphabet } : {}),
  });
}

const datum = (over: Partial<ThemeDatum> & { field?: ThemeDatum['field'] }): ThemeDatum => ({
  field: () => undefined,
  value: Number.NaN,
  letter: '',
  start: 0,
  end: 0,
  ...over,
});

const fillOf = (channel: Channel, d: ThemeDatum) => hex(channel.constant ?? channel.at(d));

function fieldReader(fields: readonly RangesField[] | undefined, row: number): ThemeDatum['field'] {
  return (name) => {
    const f = fields?.find((x) => x.name === name);
    const v = f ? cell(f.column, row) : undefined;
    return v === null ? undefined : v;
  };
}

const escape = (s: string) => s.replace(/[&<>"]/g, (c) => `&#${c.charCodeAt(0)};`);

/** Draws the section's letters, ranges, values and variants layers for a window; waits for the store's tiles. */
export async function drawSection(
  model: StaticModel,
  view: ResolvedView,
  store: DataStore,
  options: DrawOptions,
): Promise<Drawn> {
  const section: ResolvedSection | undefined = options.section
    ? view.sections.find((s) => s.address === options.section)
    : view.sections[0];
  if (!section) throw new Error('no such section');
  const axis = section.axis;
  const [ws, we] = options.window;
  const a = nearestOffsetOf(axis, ws);
  const b = nearestOffsetOf(axis, we);
  if (a === undefined || b === undefined) throw new Error('the axis waits for its length');
  const u0 = Math.min(a, b);
  const scale = (options.width - LabelWidth - 8) / (Math.abs(b - a) || 1);
  const plotX = LabelWidth;
  const plotW = options.width - LabelWidth - 8;
  // the display extent of positions `[s, e)` on the axis, a cell's edges: left, right (undefined: not shown)
  const extent = (s: number, e: number): [number, number] | undefined => {
    const first = cellOffset(axis, axis.dir === 1 ? s : e - 1);
    const last = cellOffset(axis, axis.dir === 1 ? e - 1 : s);
    return first === undefined || last === undefined ? undefined : [first, last + 1];
  };
  const xOf = (u: number) => plotX + (u - u0) * scale;

  const marks: Mark[] = [];
  const body: string[] = [];
  const skipped: string[] = [];
  let pending = 0;
  let y = Gap;

  for (const track of section.tracks) {
    const top = y;
    let height = track.height ?? 20;
    const trackMarks: Mark[] = [];
    for (const layer of track.layers) {
      const access = layerSource(model, layer);
      // a layer reached through a path (a follower's) is carried by the path's maps: left to a consumer that builds them
      if (!access || !Kinds.has(layer.type) || layer.path) {
        skipped.push(layer.address ?? track.address);
        continue;
      }
      const address = layer.address ?? track.address;
      const read = await readTables(store, access, ws, we);
      pending += read.pending;
      const add = (m: Omit<Mark, 'layer'>) => trackMarks.push({ ...m, layer: address });
      const ex = (s: number, e: number) => {
        const r = extent(s, e);
        return r ? { x: xOf(r[0]), width: (r[1] - r[0]) * scale } : undefined;
      };
      if (layer.type === 'letters') {
        const alphabet = model.systems.get(layer.data?.system ?? '')?.alphabet ?? 'protein';
        const body = { kind: 'scheme', params: { name: lettersDefault(alphabet) } } as unknown as ThemeBody;
        const color = channelOf(layer, 'color', body, alphabet);
        for (const t of read.tables) {
          const l = t.letters;
          if (!l) continue;
          const n = l.letters.length;
          for (let p = Math.max(ws, l.start); p < Math.min(we, l.start + n); p++) {
            const stored = typeof l.letters === 'string' ? l.letters.charCodeAt(p - l.start) : l.letters[p - l.start];
            // a nucleotide read against its strand on screen shows its complement (`complement`, `auto`: when the axis is reversed)
            const flip =
              alphabet !== 'protein' &&
              (layer.params.complement === 'always' || (layer.params.complement !== 'never' && axis.dir === -1));
            const c = String.fromCharCode(flip ? complementCode(stored, alphabet === 'rna') : stored);
            const at = ex(p, p + 1);
            if (at)
              add({
                kind: 'letter',
                position: p,
                end: p + 1,
                ...at,
                fill: fillOf(color, datum({ letter: c, start: p, end: p + 1 })),
                text: c,
              });
          }
        }
        height = Math.max(height, 20);
      } else if (layer.type === 'ranges' || layer.type === 'variants') {
        const color = channelOf(
          layer,
          'color',
          `token:${layer.type === 'ranges' ? 'range' : 'variant'}` as unknown as ThemeBody,
        );
        const rows: { s: number; e: number; row: number; t: Extract<Tables['ranges'], { form: 'list' }> }[] = [];
        // an item's key is its id, else its row (`_spec` Render, Data): a key a tile before showed is that range again
        const seen = new Set<string | number>();
        for (const t of read.tables) {
          const r = t.ranges;
          if (r?.form !== 'list') continue;
          for (let i = 0; i < r.start.length; i++) {
            if (r.end[i] <= ws || r.start[i] >= we) continue;
            const key = r.id?.get(i) ?? r.row?.[i];
            if (key !== undefined && key !== null) {
              if (seen.has(key)) continue;
              seen.add(key);
            }
            rows.push({ s: r.start[i], e: r.end[i], row: i, t: r });
          }
        }
        rows.sort((p, q) => p.s - q.s || p.e - q.e);
        const lanes: number[] = [];
        for (const { s, e, row, t } of rows) {
          const field = fieldReader(t.fields, row);
          const fill = fillOf(color, datum({ field, start: s, end: e }));
          const label = field('label') ?? t.id?.get(row) ?? undefined;
          const text = label === undefined ? undefined : String(label);
          const at = ex(Math.max(s, ws), Math.min(e, we));
          if (!at) continue;
          if (layer.type === 'variants') {
            add({ kind: 'variant', position: s, end: e, ...at, fill, ...(text ? { text } : {}) });
            continue;
          }
          let lane = lanes.findIndex((end) => end <= s);
          if (lane < 0) lane = lanes.push(0) - 1;
          lanes[lane] = e;
          add({ kind: 'range', position: s, end: e, ...at, fill, ...(text ? { text } : {}), value: lane });
        }
        height = layer.type === 'ranges' ? Math.max(lanes.length, 1) * LaneHeight + 6 : Math.max(height, 32);
      } else {
        const color = channelOf(layer, 'color', 'token:value' as unknown as ThemeBody);
        const vals: { p: number; v: number }[] = [];
        for (const t of read.tables) {
          const v = t.values;
          if (v?.form !== 'positions') continue;
          for (let p = Math.max(ws, v.start); p < Math.min(we, v.start + v.value.length); p++)
            vals.push({ p, v: v.value[p - v.start] });
        }
        const finite = vals.map((x) => x.v).filter((v) => !Number.isNaN(v));
        const given = layer.params.domain as readonly number[] | undefined;
        // the open domain: the values in view, down to zero
        const lo = given ? given[0] : Math.min(0, ...finite);
        const hi = given ? given[1] : Math.max(...finite);
        for (const { p, v } of vals) {
          const at = ex(p, p + 1);
          if (!at || Number.isNaN(v)) continue;
          const share = hi > lo ? Math.max(0, Math.min(1, (v - lo) / (hi - lo))) : 0;
          add({
            kind: 'value',
            position: p,
            end: p + 1,
            ...at,
            fill: fillOf(color, datum({ value: v, start: p, end: p + 1 })),
            value: v,
            share,
          });
        }
        height = Math.max(height, 32);
      }
    }
    // the track: its label, a well, its marks
    body.push(`<text x="8" y="${top + 13}" font-size="12" fill="#6b655c">${escape(track.label)}</text>`);
    body.push(`<rect x="${plotX}" y="${top}" width="${plotW}" height="${height}" fill="#e9e6e0"/>`);
    for (const m of trackMarks) body.push(markSvg(m, top, height, plotX, plotW));
    marks.push(...trackMarks);
    y = top + height + Gap;
  }
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${options.width}" height="${y}" viewBox="0 0 ${options.width} ${y}">${body.join('')}</svg>`;
  return { svg, marks, plot: { x: plotX, width: plotW }, pending, skipped };
}

const f = (n: number) => String(Math.round(n * 100) / 100);

/** A mark as SVG, clipped to the plot by its bounds. */
function markSvg(m: Mark, top: number, height: number, plotX: number, plotW: number): string {
  const x0 = Math.max(m.x, plotX);
  const x1 = Math.min(m.x + m.width, plotX + plotW);
  const w = x1 - x0;
  switch (m.kind) {
    case 'letter':
      return (
        `<rect x="${f(x0)}" y="${top}" width="${f(w)}" height="${height}" fill="${m.fill}"/>` +
        (w >= LetterCell
          ? `<text x="${f(m.x + m.width / 2)}" y="${top + height / 2 + 4}" text-anchor="middle" font-size="12" fill="#2d271e">${escape(m.text ?? '')}</text>`
          : '')
      );
    case 'range': {
      const y = top + 3 + (m.value ?? 0) * LaneHeight;
      return (
        `<rect x="${f(x0)}" y="${y}" width="${f(w)}" height="${LaneHeight - 2}" rx="3" fill="${m.fill}"/>` +
        (m.text ? `<text x="${f(x0 + 3)}" y="${y + 10}" font-size="10" fill="#2d271e">${escape(m.text)}</text>` : '')
      );
    }
    case 'value': {
      const h = (m.share ?? 0) * height;
      return `<rect x="${f(x0)}" y="${f(top + height - h)}" width="${f(w)}" height="${f(h)}" fill="${m.fill}"/>`;
    }
    case 'variant': {
      const cx = m.x + m.width / 2;
      return `<line x1="${f(cx)}" y1="${top + 8}" x2="${f(cx)}" y2="${top + height}" stroke="${m.fill}"/><circle cx="${f(cx)}" cy="${top + 8}" r="4" fill="${m.fill}"/>`;
    }
  }
}
