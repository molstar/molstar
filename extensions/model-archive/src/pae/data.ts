/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The predicted aligned error (PAE) of a model as the PAE card draws it (parity.md §6, "PAE plot"): a dense square of
// numbers over the residues the metric covers, and its picture. Plain data, no plugin and no DOM: the drawing fills an
// RGBA buffer the card puts on a canvas in one call.

/** A residue of the matrix, by author numbering (what targets name, targets.md §1). */
export interface PaeResidue {
  readonly chain: string;
  readonly seq: number;
  /** The insertion code, or ''. */
  readonly insCode: string;
  /** The residue's name (ALA). */
  readonly name: string;
}

export interface PaeMatrix {
  /** The metric's name in the file ("PAE"). */
  readonly name: string;
  /** Residues on each side. */
  readonly size: number;
  /** Row-major: `values[i * size + j]` is the error of residue j aligned on residue i (Å); NaN where the file has none. */
  readonly values: Float32Array;
  /** The highest value: the light end of the scale. */
  readonly max: number;
  readonly residues: readonly PaeResidue[];
}

/** The scale: AlphaFold DB's greens (the classic plot's), dark for a low error, near white for the highest. */
export const PaeColors = { low: [0x00, 0x44, 0x1b], high: [0xf7, 0xfc, 0xf5] } as const;
export const PaeColorStyles = { low: '#00441b', high: '#f7fcf5' } as const;

/**
 * The picture, one pixel per cell: x is residue i (the row), y residue j. A missing cell takes its mirror when that has
 * a value (a file may give one triangle); else it stays transparent. `out` is reused when it fits.
 */
export function paeImage(m: PaeMatrix, out?: Uint8ClampedArray): Uint8ClampedArray {
  const { size, values } = m;
  const length = size * size * 4;
  const rgba = out && out.length === length ? out : new Uint8ClampedArray(length);
  const scale = m.max > 0 ? 1 / m.max : 0;
  const [r0, g0, b0] = PaeColors.low;
  const [r1, g1, b1] = PaeColors.high;
  for (let i = 0; i < size; i++) {
    for (let j = 0; j < size; j++) {
      let v = values[i * size + j];
      if (Number.isNaN(v)) v = values[j * size + i];
      const p = (j * size + i) * 4;
      if (Number.isNaN(v)) {
        rgba[p + 3] = 0;
        continue;
      }
      const t = Math.min(1, Math.max(0, v * scale));
      rgba[p] = r0 + (r1 - r0) * t;
      rgba[p + 1] = g0 + (g1 - g0) * t;
      rgba[p + 2] = b0 + (b1 - b0) * t;
      rgba[p + 3] = 255;
    }
  }
  return rgba;
}

/** The cell under a point, from its position in the plot (0..1 on each axis); undefined outside. */
export function paeCellAt(m: PaeMatrix, fx: number, fy: number): { i: number; j: number } | undefined {
  if (!(fx >= 0 && fx < 1 && fy >= 0 && fy < 1)) return undefined;
  return { i: Math.floor(fx * m.size), j: Math.floor(fy * m.size) };
}

/** The value of a cell (its mirror when the file gives one triangle); NaN when there is none. */
export function paeValue(m: PaeMatrix, i: number, j: number): number {
  const v = m.values[i * m.size + j];
  return Number.isNaN(v) ? m.values[j * m.size + i] : v;
}

/** A residue as a target's address: `A:120`, `B:120A`. */
export const paeAddress = (r: PaeResidue) => `${r.chain}:${r.seq}${r.insCode}`;

/** The two residues of a cell as one target (`A:12, A:40`; one when they are the same). */
export function paePairTarget(m: PaeMatrix, i: number, j: number): string {
  const a = paeAddress(m.residues[i]);
  const b = paeAddress(m.residues[j]);
  return a === b ? a : `${a}, ${b}`;
}

/** "A 12 ALA · A 40 GLY · 3.2 Å", what the card says of the cell under the pointer. */
export function paeCellLabel(m: PaeMatrix, i: number, j: number): string {
  const a = m.residues[i];
  const b = m.residues[j];
  const v = paeValue(m, i, j);
  const value = Number.isNaN(v) ? 'no value' : `${v.toFixed(1)} Å`;
  return `${a.chain} ${a.seq}${a.insCode} ${a.name} · ${b.chain} ${b.seq}${b.insCode} ${b.name} · ${value}`;
}

/** The chains along the axis, as runs of the matrix (for the ticks): `[label, from, to)`. */
export function paeChains(m: PaeMatrix): { label: string; from: number; to: number }[] {
  const out: { label: string; from: number; to: number }[] = [];
  for (let k = 0; k < m.size; k++) {
    const chain = m.residues[k].chain;
    const last = out[out.length - 1];
    if (last && last.label === chain) last.to = k + 1;
    else out.push({ label: chain, from: k, to: k + 1 });
  }
  return out;
}
