/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ResidueSchemes } from '@molstar/core/util/color/residue-schemes';

// Where each built-in scheme's colors come from (_spec/README.md "Themes", Schemes): what a picker and the docs show.
// Kept out of the schemes themselves, so drawing (headless, an export) carries no citation it never shows.

const ColumnMeasureNote = 'needs the rows layer’s column measures';
const Lists = 'Mol*’s `ColorLists`';

const Seq: Readonly<Record<string, string>> = {
  classes: 'Seq*’s residue class tokens (`seq-aa-*`)',
  'clustal-x': `${ResidueSchemes['clustal-x'].provenance}; ${ColumnMeasureNote}`,
  hydrophobicity: 'Wimley and White 1996 (doi:10.1038/nsb1096-842), Mol*’s `ResidueHydrophobicity`, red to blue',
  'percent-identity': `Jalview’s percentage identity, blue shades by the column’s identity; ${ColumnMeasureNote}`,
  blosum62: `Jalview’s BLOSUM62 score against the column’s consensus; ${ColumnMeasureNote}`,
  espript: `ESPript (Robert and Gouet 2014): identical columns filled, similar residues by groups (an approximation of its score); ${ColumnMeasureNote}`,
  'nucleotide-classes': 'Seq*’s nucleotide tokens (`seq-nt-*`)',
  ramp: 'Seq*’s ramp tokens (`seq-ramp-*`), monotone in lightness',
  diverging: 'Seq*’s diverging tokens (`seq-div-*`), centered on 0 unless the form says `mid`',
  viridis: `Matplotlib’s viridis (van der Walt and Smith), ${Lists}`,
  cividis: `Nuñez, Anderton and Renslow 2018, ${Lists}`,
  magma: `Matplotlib’s magma, ${Lists}`,
  blues: `ColorBrewer Blues (Brewer), ${Lists}`,
  'rd-bu': `ColorBrewer RdBu (Brewer), ${Lists}`,
  plddt: 'AlphaFold DB’s model confidence bands (Jumper et al. 2021): <50, 50–70, 70–90, ≥90',
  categories: 'Seq*’s category tokens (`seq-cat-*`), checked for color-vision deficiencies',
  'okabe-ito': `Okabe and Ito 2008, ${Lists}`,
  tableau10: `Tableau 10, ${Lists}`,
  category10: `D3’s category10, ${Lists}`,
  secondary: 'Seq*’s secondary-structure tokens (`seq-ss-*`)',
  giemsa:
    'UCSC’s cytoband stains as greys (darker for stronger stain, flipped in the dark theme); `acen` in `seq-centromere`',
};

/** A built-in scheme's provenance (a residue scheme's is core's, with its data); undefined for one Seq* did not build in. */
export function schemeProvenance(name: string): string | undefined {
  return (
    Seq[name] ??
    (Object.hasOwn(ResidueSchemes, name) ? ResidueSchemes[name as keyof typeof ResidueSchemes].provenance : undefined)
  );
}
