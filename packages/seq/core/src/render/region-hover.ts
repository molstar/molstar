/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a hover over data read by region says (_spec/README.md "Render", Hover texts by region): a pairs bin as its two
// bin ranges, its count and its level, never as a pair; a region grid's cell as its bin on each axis and its value (a
// sparse tile's unlisted cell its `fill`, a NaN "no data", a tile not loaded "loading"); a matrix cell as its item, its
// state and its value or genotype, at a coarse level its bin and its mean or majority with its share. Texts only: a
// host's card lays them out (the viewer's hover card), nothing is drawn. Read on a hover, from what the frame drew.

import { spanText } from './accessibility.js';
import type { GridLay } from './grid.js';
import { regionCellAt } from './grid-region.js';
import { columnsOf, type MatrixLay } from './matrix.js';
import { aOf, bOf, type PairsLay } from './pairs-lay.js';

/** A hover's words: a title, a card's rows (label, value), and the sentence an announcer says. */
export interface RegionHover {
  readonly title: string;
  readonly rows: readonly (readonly [string, string])[];
  readonly text: string;
}

/** A system's name as a card says it (the viewer's: a sequence's name, a chain's label). */
export type SystemName = (system: string) => string;

/** A bin's width in positions as a card says it: `256 bp bins`, `5 kb bins`, `1 Mb bins`. */
export function binsText(level: number): string {
  const unit = (v: number, u: string) => `${Number(v.toFixed(v < 10 ? 1 : 0))} ${u} bins`;
  if (level >= 1e6) return unit(level / 1e6, 'Mb');
  if (level >= 1e3) return unit(level / 1e3, 'kb');
  return `${level.toLocaleString('en-US')} bp bins`;
}

/** A genotype class as VCF writes it (`0/0`, `0/1`, `1/1`, a partly called `./1`, an uncalled `./.`). */
export function genotypeText(code: number): string {
  switch (code) {
    case 0:
      return '0/0';
    case 1:
      return '0/1';
    case 2:
      return '1/1';
    case 3:
      return './1';
    default:
      return './.';
  }
}

function numberText(v: number): string {
  if (!Number.isFinite(v)) return String(v);
  if (Number.isInteger(v)) return v.toLocaleString('en-US');
  const a = Math.abs(v);
  return Number(v.toFixed(a >= 100 ? 1 : a >= 1 ? 2 : 3)).toString();
}

function hover(title: string, rows: (readonly [string, string])[], words: string[]): RegionHover {
  return { title, rows, text: words.join(', ') };
}

/**
 * A pairs layer's bin pair by its element key (`bin:…`, `pairs-region.ts`): its two bin ranges, its count, its level;
 * undefined for a key that is not a bin drawn.
 */
export function pairsBinHover(lay: PairsLay, key: string | number, name: SystemName): RegionHover | undefined {
  const prep = lay.prep;
  const d = prep?.data;
  if (!prep || !d || prep.level === undefined || typeof key !== 'string' || !key.startsWith('bin:')) return undefined;
  const k = prep.index.get(key);
  if (k === undefined) return undefined;
  const a = `${name(aOf(prep, k))}:${spanText(d.aStart[k], d.aEnd[k])}`;
  const b = `${name(bOf(prep, k))}:${spanText(d.bStart[k], d.bEnd[k])}`;
  const n = prep.weight?.[k] ?? 1;
  const count = `${numberText(n)} pair${n === 1 ? '' : 's'}`;
  const bins = binsText(prep.level);
  const reversed = d.bStrand[k] < 0 || d.aStrand[k] < 0;
  const kind = d.kinds[d.kind[k]];
  const rows: (readonly [string, string])[] = [
    ['From', a],
    ['To', b],
    ['Pairs', kind ? `${count} (${kind})` : count],
    ['Bins', bins],
  ];
  if (reversed) rows.push(['Strand', 'reversed']);
  return hover(`${a} ↔ ${b}`, rows, [
    `${a} to ${b}`,
    count,
    ...(kind ? [kind] : []),
    ...(reversed ? ['reversed'] : []),
    bins,
    'zoom in to see them',
  ]);
}

/**
 * A region grid's cell under (i, j): its bin on each axis and its value (`count` and such as the layer's field), a
 * sparse tile's unlisted cell its `fill`, NaN "no data", a tile not loaded "loading"; undefined for a held grid.
 */
export function gridCellHover(lay: GridLay, i: number, j: number, name: SystemName): RegionHover | undefined {
  if (!lay.region) return undefined;
  const L = lay.binLevel;
  const x0 = Math.floor(i / L) * L;
  const y0 = Math.floor(j / L) * L;
  const x = `${name(lay.system)}:${spanText(x0, x0 + L)}`;
  const y = `${name(lay.ySystem)}:${spanText(y0, y0 + L)}`;
  const v = regionCellAt(lay, i, j);
  const value = v === undefined ? 'loading' : Number.isNaN(v) ? 'no data' : numberText(v);
  const rows: (readonly [string, string])[] = [
    ['x', x],
    ['y', y],
    ['Value', value],
  ];
  if (L > 1) rows.push(['Bins', binsText(L)]);
  return hover(`${x} × ${y}`, rows, [`${x} × ${y}`, `value ${value}`, ...(L > 1 ? [binsText(L)] : [])]);
}

/**
 * A matrix cell under position `i` and state `s`: its item (its id and extent, the shortest item there) and its value
 * or genotype; at a coarse level its bin and its mean, or its majority class with its share; undefined off the data.
 */
export function matrixCellHover(lay: MatrixLay, i: number, s: number, name: SystemName): RegionHover | undefined {
  if (s < 0 || s >= lay.states) return undefined;
  const state = lay.names[s] ?? `state ${(s + 1).toLocaleString('en-US')}`;
  const system = name(lay.system);
  for (const t of lay.tables) {
    if (t.form === 'items') {
      let best = -1;
      for (let k = 0; k < t.start.length; k++) {
        const e = Math.max(t.end[k], t.start[k] + 1);
        if (t.start[k] <= i && i < e && (best < 0 || e - t.start[k] < t.end[best] - t.start[best])) best = k;
      }
      if (best < 0) continue;
      const at = `${system}:${spanText(t.start[best], Math.max(t.end[best], t.start[best] + 1))}`;
      const id = t.id?.get(best);
      const item = id ? `${id} (${at})` : at;
      const c = best * t.states + s;
      const value = t.codes ? genotypeText(t.codes[c]) : t.value ? valueOf(t.value[c]) : '';
      return hover(
        item,
        [['Item', item], ['State', state], ...(value ? [['Value', value] as const] : [])],
        [item, state, ...(value ? [value] : [])],
      );
    }
    const n = columnsOf(t);
    const level = lay.level;
    const b = Math.floor((i - t.start) / level);
    if (b < 0 || b >= n) continue;
    const start = t.start + b * level;
    const at = `${system}:${spanText(start, start + level)}`;
    const c = b * t.states + s;
    const value = t.classes
      ? `${genotypeText(t.classes[c])}${t.shares ? ` (${Math.round((t.shares[c] / 255) * 100)}%)` : ''}`
      : t.mean
        ? `mean ${valueOf(t.mean[c])}`
        : '';
    return hover(
      at,
      [['Bin', at], ['State', state], ...(value ? [['Value', value] as const] : []), ['Bins', binsText(level)]],
      [at, state, ...(value ? [value] : []), binsText(level)],
    );
  }
  return undefined;
}

function valueOf(v: number): string {
  return Number.isNaN(v) ? 'no data' : numberText(v);
}
