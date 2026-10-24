/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { vmdColor } from '../colors.js';
import { fail, nearest, type Span, warning } from '../diagnostics.js';
import type { ColorSpec, RepType, VmdRep } from '../ir.js';
import type { Compiler } from '../table.js';

// VMD's drawing styles, coloring methods and materials as the IR's (vmd.md §2–3). Names are case-insensitive, as in
// VMD; a style or method Mol* does not have fails, naming the ones it does.

const Styles: Readonly<Record<string, RepType>> = {
  lines: 'lines',
  bonds: 'sticks',
  licorice: 'sticks',
  cpk: 'ball-and-stick',
  vdw: 'spheres',
  points: 'points',
  newcartoon: 'cartoon',
  cartoon: 'cartoon',
  newribbons: 'cartoon',
  ribbons: 'cartoon',
  tube: 'tube',
  trace: 'backbone',
  quicksurf: 'gaussian-surface',
  msms: 'surface',
  surf: 'surface',
  newsurf: 'surface',
};

/** The styles' names as VMD writes them, for help and completion. */
export const VmdStyleNames = [
  'Lines',
  'Bonds',
  'Licorice',
  'CPK',
  'VDW',
  'Points',
  'NewCartoon',
  'Cartoon',
  'NewRibbons',
  'Ribbons',
  'Tube',
  'Trace',
  'QuickSurf',
  'MSMS',
  'Surf',
  'NewSurf',
];

const UnsupportedStyles: Readonly<Record<string, string>> = {
  hbonds: 'Hydrogen bonds are shown from the Structure panel’s interactions.',
  dynamicbonds: 'Bonds come from the file and Mol*’s rules.',
  paperchain: 'Not drawn by Mol*.',
  twister: 'Not drawn by Mol*.',
  beads: 'Not drawn by Mol*.',
  dotted: 'Not drawn by Mol*: use Points.',
  solvent: 'Not drawn by Mol*: use VDW.',
  polyhedra: 'Not drawn by Mol*.',
  isosurface: 'Volumes are loaded from the Volume panel.',
  volumeslice: 'Volumes are loaded from the Volume panel.',
  fieldlines: 'Not drawn by Mol*.',
  orbital: 'Not drawn by Mol*.',
};

export function vmdStyle(text: string, at: Span, c: Compiler): RepType {
  const t = text.trim().toLowerCase();
  const style = Styles[t];
  if (style) {
    if (t === 'ribbons' || t === 'newribbons') {
      c.report(warning('unsupported-argument', `${text} shows as Mol*’s cartoon.`, at));
    }
    return style;
  }
  const why = UnsupportedStyles[t];
  if (why) fail('unsupported', `The style ${text} is not supported.`, at, why);
  const near = nearest(t, Object.keys(Styles));
  return fail(
    'argument',
    `'${text}' is not a VMD drawing style.`,
    at,
    near
      ? `Did you mean ${VmdStyleNames.find((n) => n.toLowerCase() === near)}?`
      : `Styles: ${VmdStyleNames.join(', ')}.`,
  );
}

const Methods: Readonly<Record<string, ColorSpec>> = {
  name: { theme: 'element-symbol', params: { palette: 'vmd' } },
  element: { theme: 'element-symbol', params: { palette: 'vmd' } },
  type: { theme: 'element-symbol', params: { palette: 'vmd' } },
  resname: { theme: 'residue-name' },
  resid: { theme: 'sequence-id' },
  chain: { theme: 'chain-id' },
  segname: { theme: 'chain-id' },
  structure: { theme: 'secondary-structure' },
  beta: { theme: 'uncertainty' },
  occupancy: { theme: 'occupancy' },
  charge: { theme: 'partial-charge' },
};

export const VmdColorMethodNames = [
  'ColorID',
  'Name',
  'Element',
  'Type',
  'ResName',
  'ResID',
  'Chain',
  'SegName',
  'Molecule',
  'Structure',
  'Beta',
  'Occupancy',
  'Charge',
];

/** A coloring method and its arguments (`ColorID 1`, `Chain`, `Molecule` with the molecule's id). */
export function vmdColorMethod(words: readonly string[], at: Span, molecule: number | undefined): ColorSpec {
  const t = (words[0] ?? '').trim().toLowerCase();
  if (t === 'colorid') {
    const id = words[1];
    const rgb = id === undefined ? undefined : vmdColor(id);
    if (rgb === undefined) fail('argument', `ColorID takes a color id from 0 to 32, not '${id ?? ''}'.`, at);
    return { rgb };
  }
  if (t === 'molecule') {
    // VMD's Molecule coloring: the molecule's id as a ColorID.
    return { rgb: vmdColor(String((molecule ?? 0) % 33))! };
  }
  const m = Methods[t];
  if (m) return m;
  const near = nearest(t, Object.keys(Methods));
  return fail(
    'unsupported',
    `The coloring method '${words[0] ?? ''}' is not supported.`,
    at,
    near
      ? `Did you mean ${VmdColorMethodNames.find((n) => n.toLowerCase() === near)}?`
      : `Methods: ${VmdColorMethodNames.join(', ')}.`,
  );
}

/** A material: Opaque and Transparent; the others draw as Opaque, with a warning. */
export function vmdMaterial(text: string, at: Span, c: Compiler): VmdRep['material'] {
  const t = text.trim().toLowerCase();
  if (t === 'transparent' || t === 'glass1' || t === 'glass2' || t === 'glass3' || t === 'ghost') return 'transparent';
  if (t !== 'opaque') c.report(warning('unsupported-argument', `The material ${text} draws as Opaque.`, at));
  return 'opaque';
}
