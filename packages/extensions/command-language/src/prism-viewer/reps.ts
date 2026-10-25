/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import { ElementSymbolColors } from '@molstar/graphics/theme/color/element-symbol';
import { VmdNameColors } from '../language/colors.js';
import type { ColorSpec, RepType } from '../language/ir.js';

// The IR's representations and colors as Mol*'s (pymol.md §2, vmd.md §2–3): a representation type and the params
// that make it look like the program's (PyMOL's sticks are licorice: bonds as thick as atoms), and color themes with
// the params the language asks for.

export interface MolstarRep {
  readonly type: string;
  /** Params of the type, set after the component is made (`representation.update`). */
  readonly typeParams?: Readonly<Record<string, unknown>>;
  readonly sizeTheme?: { readonly name: string; readonly params?: Readonly<Record<string, unknown>> };
}

export const MolstarReps: Readonly<Record<RepType, MolstarRep>> = {
  cartoon: { type: 'cartoon' },
  putty: { type: 'putty' },
  tube: { type: 'putty', sizeTheme: { name: 'uniform', params: { value: 1 } } },
  backbone: { type: 'backbone' },
  sticks: { type: 'ball-and-stick', typeParams: { sizeAspectRatio: 1 } },
  'ball-and-stick': { type: 'ball-and-stick' },
  lines: { type: 'line' },
  spheres: { type: 'spacefill' },
  surface: { type: 'molecular-surface' },
  mesh: { type: 'molecular-surface', typeParams: { visuals: ['molecular-surface-wireframe'] } },
  'gaussian-surface': { type: 'gaussian-surface' },
  points: { type: 'point' },
  ellipsoids: { type: 'ellipsoid' },
};

/** What a representation the console did not make shows, in the IR's words; undefined: nothing PyMOL shows. */
export function repTypeOf(
  molstarType: string,
  typeParams?: { readonly visuals?: readonly string[] },
): RepType | undefined {
  switch (molstarType) {
    case 'cartoon':
      return 'cartoon';
    case 'putty':
      return 'putty';
    case 'backbone':
      return 'backbone';
    case 'ball-and-stick':
      return 'sticks';
    case 'line':
      return 'lines';
    case 'spacefill':
      return 'spheres';
    case 'molecular-surface':
      return typeParams?.visuals?.every((v) => v.includes('wireframe')) ? 'mesh' : 'surface';
    case 'gaussian-surface':
      return 'gaussian-surface';
    case 'point':
      return 'points';
    case 'ellipsoid':
      return 'ellipsoids';
    default:
      return undefined;
  }
}

export const hex = (rgb: number) => `#${rgb.toString(16).padStart(6, '0')}`;

/** Mol*'s color list for a PyMOL palette (`spectrum`). */
const Palettes: Readonly<Record<string, readonly number[]>> = {
  rainbow: [0x0000ff, 0x00ffff, 0x00ff00, 0xffff00, 0xff0000],
  rainbow_rev: [0xff0000, 0xffff00, 0x00ff00, 0x00ffff, 0x0000ff],
  blue_white_red: [0x0000ff, 0xffffff, 0xff0000],
  red_white_blue: [0xff0000, 0xffffff, 0x0000ff],
  blue_red: [0x0000ff, 0xff0000],
  red_blue: [0xff0000, 0x0000ff],
};

/** A theme as `representation.update` takes it: its name and params. */
export function colorTheme(c: Extract<ColorSpec, { theme: string }>): {
  name: string;
  params: Record<string, unknown>;
} {
  const p = c.params ?? {};
  if (c.theme === 'element-symbol') {
    const params: Record<string, unknown> = {};
    if (p.palette === 'vmd') {
      params.carbonColor = { name: 'uniform', params: { value: Color(VmdNameColors.C) } };
      params.colors = {
        name: 'custom',
        params: Object.fromEntries(
          Object.keys(ElementSymbolColors).map((e) => [
            e,
            Color(VmdNameColors[e] ?? (ElementSymbolColors as Record<string, number>)[e]),
          ]),
        ),
      };
      params.lightness = 0;
    } else if (typeof p.carbon === 'number') {
      params.carbonColor = { name: 'uniform', params: { value: Color(p.carbon) } };
    }
    return { name: c.theme, params };
  }
  if (typeof p.palette === 'string' && Palettes[p.palette]) {
    return {
      name: c.theme,
      params: { list: { kind: 'interpolate', colors: Palettes[p.palette].map((x) => Color(x)) } },
    };
  }
  return { name: c.theme, params: {} };
}
