/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A residue's text in its cell, part of the renderer API (seqspec:viewer#the-render-backend): the one-letter font size
// and whether a row shows text at all, shared by the letters, the rows and the emphasis drawn over them, so a row never
// shows some columns' text and not others' (VL-06). Doubles are read from `letterRoom`, never passed (a double passed
// to a call is boxed below TurboFan).

import { Molecules } from '@molstar/seq-core/model';
import { AxisGeometry } from '../axis.js';
import { letterCap } from '../layout.js';
import { Lod } from '../viewport.js';
import type { DrawView } from './view.js';

// Read once: a module member read per frame allocates under the test runner's module wrapper.
const WRAPPED = AxisGeometry.wrapped;
const LETTERS = Lod.letters;

/** The row's height (`[0]`) and the font's widest letter (`[1]`; 0: not measured), CSS px, for the two below. */
export const letterRoom = new Float64Array(2);

/**
 * Block b's one-letter font px in a row `letterRoom[0]` tall: the cell's 1.25, the row's 0.72, at most the cap
 * (`letterCap`), at least 6; whole px. Every residue text of a row is drawn at it.
 */
export function letterFontPx(v: DrawView, b: number): number {
  return (
    Math.max(6, Math.min(letterCap(v.props), Math.floor(Math.min(v.cell * v.span[b] * 1.25, letterRoom[0] * 0.72)))) | 0
  );
}

/**
 * Whether block b's residues show as text in a row `letterRoom[0]` tall whose font's widest letter is `letterRoom[1]`:
 * a wrapped line always; linear at the letters LOD, a row of 10 px or more, and a cell holding that letter and
 * `minLetterCell`. Otherwise their cells (or marks) show, and nothing drawn over them shows text either.
 */
export function textShown(v: DrawView, b: number): boolean {
  if (v.geometry === WRAPPED) return true;
  const w = v.cell * v.span[b];
  return v.lod === LETTERS && letterRoom[0] >= 10 && w >= letterRoom[1] && w >= v.props.minLetterCell;
}

const SYMBOLS: string[] = [];
/** An atom's element symbol by atomic number ("X" when none), made on first use: the letters' and emphasis's. */
export function elementSymbol(z: number): string {
  SYMBOLS[z] ??= Molecules.elementSymbol(z) || 'X';
  return SYMBOLS[z];
}
