/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { DataStore, Want } from '@molstar/seqx-core/data';
import { type Location, location } from '@molstar/seqx-core/query';
import { gridSummary } from '@molstar/seqx-core/render';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedView } from '@molstar/seqx-core/view';

// The readout of a grid over a complex, as a host's panel writes it (the panel is the application's: the story draws
// one): a hovered cell named by both residues with their chains and its value; a selected rect named by its parts,
// aligned (x) and scored (y), and the PAE's mean over it (`gridSummary`), AlphaFold DB's wording.

/** What the readout reads: the view's model and resolved view (a layout's input has both), and its store. */
export interface ReadoutView {
  readonly model: StaticModel;
  readonly view: ResolvedView;
  /** The store the letters read from a source are in (none: those letters are left out). */
  readonly data?: Pick<DataStore, 'keysOf' | 'get'>;
}

/** A system's chain label: the label of the first block on it in a section of blocks, else its address. */
function chainOf(v: ReadoutView, system: string): string {
  for (const s of v.view.sections) for (const b of s.blocks ?? []) if (b.axis.system === system) return b.label;
  return system;
}

/**
 * What the readout wants of the store: the letters of each block's system read from a source, whole (a host's panel
 * holds them: a layout zoomed out reads none).
 */
export function readoutWants(v: ReadoutView): Want[] {
  const out: Want[] = [];
  const seen = new Set<string>();
  for (const sec of v.view.sections)
    for (const b of sec.blocks ?? []) {
      const s = v.model.systems.get(b.axis.system);
      if (!s?.source || s.length === undefined || seen.has(s.address)) continue;
      seen.add(s.address);
      out.push({ source: s.source, system: s.name, start: 0, end: s.length, level: 1 });
    }
  return out;
}

/** A residue's letter: inline, else in the loaded tile of its source that holds it ('' while none does). */
function letterOf(v: ReadoutView, system: string, position: number): string {
  const s = v.model.systems.get(system);
  if (s?.letters !== undefined) return s.letters[position] ?? '';
  if (!s?.source || !v.data) return '';
  for (const k of v.data.keysOf({ source: s.source, system: s.name, start: position, end: position + 1, level: 1 })) {
    const e = v.data.get(k);
    const t = e?.status === 'loaded' ? e.tables.letters : undefined;
    if (!t || position < t.start || position >= t.start + t.letters.length) continue;
    const i = position - t.start;
    // letters as text, or as codes of the table's alphabet
    return typeof t.letters === 'string' ? t.letters[i] : (t.alphabet?.[t.letters[i]] ?? '');
  }
  return '';
}

/** A residue as text: its letter and its number (from 1), `K33`. */
function residue(v: ReadoutView, system: string, position: number): string {
  return `${letterOf(v, system, position)}${position + 1}`;
}

/** The value of the grids at a cell (the first grid on its two systems holding one), NaN when none. */
function valueAt(v: ReadoutView, system: string, position: number, ySystem: string, yPosition: number): number {
  const cell = location({
    rects: [
      { a: { system, ranges: [position, position + 1] }, b: { system: ySystem, ranges: [yPosition, yPosition + 1] } },
    ],
  });
  return gridSummary(v, cell).mean;
}

/** A hovered cell: `CDK2 K33 × cyclin A2 L263: 4 Å` (scored residue first, aligned on the second: x, then y). */
export function cellReadout(
  v: ReadoutView,
  cell: { readonly system: string; readonly position: number; readonly ySystem: string; readonly yPosition: number },
): string {
  const value = valueAt(v, cell.system, cell.position, cell.ySystem, cell.yPosition);
  const a = `${chainOf(v, cell.system)} ${residue(v, cell.system, cell.position)}`;
  const b = `${chainOf(v, cell.ySystem)} ${residue(v, cell.ySystem, cell.yPosition)}`;
  return `${a} × ${b}${Number.isNaN(value) ? '' : `: ${Math.round(value * 10) / 10} Å`}`;
}

/** Each part's chain and range, 1-based: `CDK2 1–298, cyclin A2 173–432`. */
function partsText(v: ReadoutView, parts: readonly { readonly system: string; readonly ranges: ArrayLike<number> }[]) {
  const out: string[] = [];
  for (const p of parts)
    for (let i = 0; i + 1 < p.ranges.length; i += 2)
      out.push(`${chainOf(v, p.system)} ${p.ranges[i] + 1}–${p.ranges[i + 1]}`);
  return out.join(', ');
}

/**
 * A selected rect: `Aligned CDK2 1–298, Scored cyclin A2 173–432, mean 10.1 Å` (AlphaFold DB's wording: x the
 * residues aligned on, y those scored), each side's parts once; empty when the selection holds no rect.
 */
export function rectReadout(v: ReadoutView, loc: Location): string {
  if (!loc.rects.length) return '';
  const xs = new Map<string, { system: string; ranges: number[] }>();
  const ys = new Map<string, { system: string; ranges: number[] }>();
  for (const r of loc.rects) {
    const add = (m: typeof xs, p: (typeof r)['a']) => {
      const key = `${p.system}|${Array.from(p.ranges).join(',')}`;
      if (!m.has(key)) m.set(key, { system: p.system, ranges: Array.from(p.ranges) });
    };
    add(xs, r.a);
    add(ys, r.b);
  }
  const s = gridSummary(v, loc);
  const mean = Number.isNaN(s.mean) ? '' : `, mean ${s.mean.toFixed(1)} Å`;
  return `Aligned ${partsText(v, [...xs.values()])}, Scored ${partsText(v, [...ys.values()])}${mean}`;
}
