/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { channelLegends, type Legend, type PairsSource } from '@molstar/seqx-core/extension';
import { ColorMode } from './colour.js';
import type { RondoLay } from './lay.js';
import { CardLines, links, valueText } from './words.js';

// The card's lines on kinds of links (_spec/README.md, Kinds of links, The card): a chunk `rondoViewer` loads once
// composed (as the stem menu's), installed into `CardLines`: the headless entry and a first viewer chunk pay nothing.

/**
 * The fields a ribbon's links share (_spec/README.md, Kinds of links, The card): their kind, their ends' strands, then
 * every other field the data carries (`q`, gene names where the file has them): a value they all share, a number's
 * range where they differ (`q: 0.001–0.04`), absent values and differing text left out. One link: its own fields.
 */
export function fieldLines(lay: RondoLay, r: number): string[] {
  const read = lay.read as PairsSource | undefined;
  const d = read?.data;
  const datum = read?.datum;
  if (!d || !datum) return [];
  const rows: number[] = [];
  for (let k = 0; k < d.n; k++) if (lay.rowRib[k] === r) rows.push(k);
  if (!rows.length) return [];
  const shared = (name: string): string | number | undefined => {
    let v: string | number | undefined;
    let lo = Number.POSITIVE_INFINITY;
    let hi = Number.NEGATIVE_INFINITY;
    let same = true;
    for (let i = 0; i < rows.length; i++) {
      const x = datum.set(rows[i]).field(name);
      if (x === undefined || x === null || x === '' || (typeof x === 'number' && Number.isNaN(x))) return undefined;
      if (i === 0) v = x;
      else if (x !== v) same = false;
      if (typeof x === 'number') {
        lo = Math.min(lo, x);
        hi = Math.max(hi, x);
      }
    }
    if (same) return v;
    return typeof v === 'number' ? `${valueText(lo)}–${valueText(hi)}` : undefined;
  };
  const out: string[] = [];
  const type = shared('type');
  if (type !== undefined) out.push(`Type: ${type}`);
  const sa = shared('a_strand');
  const sb = shared('b_strand');
  // unstranded ends (`.`, a contact's) say none
  const sign = (x: string | number) => (x === '-' ? '−' : x === '+' ? '+' : '.');
  if (sa !== undefined && sb !== undefined && (sa !== '.' || sb !== '.')) out.push(`Strands: ${sign(sa)} ${sign(sb)}`);
  for (const name of Object.keys(d.extra ?? {})) {
    if (name === 'count') continue;
    const v = shared(name);
    if (v !== undefined) out.push(`${name}: ${valueText(v)}`);
  }
  return out;
}

/**
 * Each layer's links on the section between two ranges (`b` absent: with an end in `a`), one line a layer
 * (`Contacts: 140 links`) when the section has more than one laid out: a bin's or a stem's count per kind. Read from
 * each layer's own stems and ribbons as laid out (their links counted at layout), those meeting the ranges: bounded by
 * what the layers drew, no row read again, nothing made but the lines.
 */
export function perLayer(lay: RondoLay, sysA: string, a0: number, a1: number, sysB?: string, b0 = 0, b1 = 0): string[] {
  if (lay.peers.length < 2) return [];
  const out: string[] = [];
  for (let i = 0; i < lay.peers.length; i++) {
    const p = lay.peers[i];
    if (!p.read || (p.empty && p !== lay)) continue;
    let n = 0;
    if (sysB === undefined) {
      for (let s = 0; s < p.ns; s++)
        if (p.systems[p.stemSystem[s]] === sysA && p.stemStart[s] < a1 && p.stemEnd[s] > a0)
          n += p.stemLinks[s] + p.withinLinks[s];
    } else
      for (let r = 0; r < p.nr; r++) {
        const x = p.ribA[r];
        const y = p.ribB[r];
        const xa = p.systems[p.stemSystem[x]] === sysA && p.stemStart[x] < a1 && p.stemEnd[x] > a0;
        const yb = p.systems[p.stemSystem[y]] === sysB && p.stemStart[y] < b1 && p.stemEnd[y] > b0;
        const xb = p.systems[p.stemSystem[x]] === sysB && p.stemStart[x] < b1 && p.stemEnd[x] > b0;
        const ya = p.systems[p.stemSystem[y]] === sysA && p.stemStart[y] < a1 && p.stemEnd[y] > a0;
        if ((xa && yb) || (xb && ya)) n += p.ribLinks[r];
      }
    out.push(`${p.label || 'Links'}: ${links(n)}`);
  }
  return out;
}

/**
 * The layer's legend (_spec/README.md, Colour, The legend): a field colour its channel's values; a fixed colour one
 * swatch named by the layer (its `label`, else its track's); a named layer in the blocks' hues one swatch saying so;
 * an unnamed one in the hues none (the hue ring keys them).
 */
export function legendRondo(lay: RondoLay, options: { readonly title: string }): readonly Legend[] {
  // its widths' scale: a CSS px a measure (links, or the field summed), its own or the section's
  const per = lay.k > 0 ? Number((1 / lay.k).toPrecision(2)) : 0;
  const scale =
    per > 0
      ? `1 px = ${per.toLocaleString('en-US')} ${lay.width === 'count' ? (per === 1 ? 'link' : 'links') : lay.width}${lay.shared ? ', shared' : ''}`
      : '';
  // a field's values, unless these bins carry none (drawn in one colour): then the legend says that
  if (lay.mode === ColorMode.field && !lay.inkBins)
    return channelLegends([lay.channel], {
      title: scale ? `${options.title} (${scale})` : options.title,
      layer: lay.key,
    });
  if (lay.mode === ColorMode.system && !lay.label) return [];
  const fixed = lay.mode !== ColorMode.system;
  const name = lay.label || options.title;
  const what = lay.inkBins
    ? `${name} (bins in one colour: ${lay.channel?.field ?? 'the colour'} read link by link)`
    : fixed
      ? name
      : `${name}: by chromosome`;
  const label = scale ? `${what} · ${scale}` : what;
  // by chromosome: the words alone (a transparent swatch), never one hue for a layer drawn in all of them
  const color = fixed ? lay.one : 0;
  return [
    {
      key: `rondo:${lay.key}`,
      title: options.title,
      channel: 'color',
      entries: { kind: 'swatches', items: [{ label, color: fixed ? alpha(color, 1) : 0 }] },
      open: false,
      source: 'vocabulary',
      more: 0,
      notes: [],
      layers: [lay.key],
    },
  ];
}

/** `0xRRGGBBAA` at an opacity. */
const alpha = (c: number, a: number) => ((c & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;

/** Puts the lines into the describe hook's reach (`CardLines`). */
export function install(): void {
  CardLines.fields = fieldLines;
  CardLines.perLayer = perLayer;
  CardLines.legend = legendRondo;
}
