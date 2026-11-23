/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PairsSource } from '@molstar/seqx-core/extension';
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
  if (sa !== undefined && sb !== undefined) out.push(`Strands: ${sa === '-' ? '−' : '+'} ${sb === '-' ? '−' : '+'}`);
  for (const name of Object.keys(d.extra ?? {})) {
    if (name === 'count') continue;
    const v = shared(name);
    if (v !== undefined) out.push(`${name}: ${valueText(v)}`);
  }
  return out;
}

/**
 * Each layer's links on the section between two ranges (`b` absent: with an end in `a`), one line a layer
 * (`Contacts: 140 links`) when the section has more than one: a bin's or a stem's count per kind. Counted from each
 * layer's rows as read (bins' starts, items' ends), each row once.
 */
export function perLayer(lay: RondoLay, sysA: string, a0: number, a1: number, sysB?: string, b0 = 0, b1 = 0): string[] {
  if (lay.peers.length < 2) return [];
  const out: string[] = [];
  for (const p of lay.peers) {
    const read = p.read;
    const d = read?.data;
    if (!read || !d || p.empty) continue;
    const sysOf = (unit: Int32Array | undefined, k: number, side: 0 | 1) =>
      unit ? p.systems[unit[k]] : p.systems[Math.min(side, p.systems.length - 1)];
    let n = 0;
    for (let k = 0; k < d.n; k++) {
      if (!read.shown[k]) continue;
      const sa = sysOf(d.aUnit, k, 0);
      const sb = sysOf(d.bUnit, k, 1);
      const inA = (s: string, x: number) => s === sysA && x >= a0 && x < a1;
      const inB = (s: string, x: number) => s === sysB && x >= b0 && x < b1;
      const hit =
        sysB === undefined
          ? inA(sa, d.aStart[k]) || inA(sb, d.bStart[k])
          : (inA(sa, d.aStart[k]) && inB(sb, d.bStart[k])) || (inA(sb, d.bStart[k]) && inB(sa, d.aStart[k]));
      if (hit) n += read.weight ? read.weight[k] : 1;
    }
    out.push(`${p.label || 'Links'}: ${links(n)}`);
  }
  return out;
}

/** Puts the lines into the describe hook's reach (`CardLines`). */
export function install(): void {
  CardLines.fields = fieldLines;
  CardLines.perLayer = perLayer;
}
