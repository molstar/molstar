/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  isSaccharideShapeDivided,
  SaccharideColors,
  SaccharideCompIdMap,
} from '@molstar/model/model/structure/structure/carbohydrates/constants';
import { ElementSymbolColors } from '@molstar/graphics/theme/color/element-symbol';
import { ResidueNameColors } from '@molstar/graphics/theme/color/residue-name';
import { SecondaryStructureColors as SecStrColors } from '@molstar/graphics/theme/color/secondary-structure';
import { Color } from '@molstar/core/util/color';
import { ValuePaletteLists } from '@molstar/core/util/color/value-palette';
import { mapObjectMap, omitObjectKeys } from '@molstar/core/util/object';
import { ColorDictNameT, ColorListNameT, ColorT } from '@molstar/mvs-builder/tree/mvs/param-types';

/** Colors for individual elements, except carbon */
const ElementSymbolColorDict: Record<string, ColorT> = MvsColorDict(omitObjectKeys(ElementSymbolColors, ['C'])); // ommitting carbon color to allow easier combination of multiple color layers

/** Colors for individual amino acids, based on Jmol */
const ResidueNameColorDict: Record<string, ColorT> = MvsColorDict(ResidueNameColors);

/** Colors for amino acid groups, based on Clustal (https://www.jalview.org/help/html/colourSchemes/clustal.html) */
const AminoGroupColors: Record<string, ColorT> = {
  aromatic: '#15A4A4',
  hydrophobic: '#80A0F0',
  polar: '#15C015',
  positive: '#F01505',
  negative: '#C048C0',
  proline: '#C0C000',
  cysteine: '#F08080',
  glycine: '#F09048',
};

/** Colors for individual amino acids, based on Clustal (https://www.jalview.org/help/html/colourSchemes/clustal.html), plus Jmol colors for nucleotides (http://jmol.sourceforge.net/jscolors/) */
const ResiduePropertyColorDict: Record<string, ColorT> = {
  ...ResidueNameColorDict, // to keep nucleotide coloring
  HIS: AminoGroupColors.aromatic,
  TYR: AminoGroupColors.aromatic,
  ALA: AminoGroupColors.hydrophobic,
  VAL: AminoGroupColors.hydrophobic,
  LEU: AminoGroupColors.hydrophobic,
  ILE: AminoGroupColors.hydrophobic,
  MET: AminoGroupColors.hydrophobic,
  PHE: AminoGroupColors.hydrophobic,
  TRP: AminoGroupColors.hydrophobic,
  SER: AminoGroupColors.polar,
  THR: AminoGroupColors.polar,
  ASN: AminoGroupColors.polar,
  GLN: AminoGroupColors.polar,
  LYS: AminoGroupColors.positive,
  ARG: AminoGroupColors.positive,
  ASP: AminoGroupColors.negative,
  GLU: AminoGroupColors.negative,
  PRO: AminoGroupColors.proline,
  CYS: AminoGroupColors.cysteine,
  GLY: AminoGroupColors.glycine,
};

/** Colors for secondary structure types, based on Jmol colors (http://jmol.sourceforge.net/jscolors/) */
const SecondaryStructureColors: Record<string, ColorT> = MvsColorDict({
  // Simple categories
  helix: SecStrColors.alphaHelix,
  strand: SecStrColors.betaStrand,
  turn: SecStrColors.betaTurn,
  bend: SecStrColors.bend,

  // DSSP categories
  H: SecStrColors.alphaHelix,
  B: SecStrColors.betaStrand,
  E: SecStrColors.betaStrand,
  G: SecStrColors.threeTenHelix,
  I: SecStrColors.piHelix,
  P: Color(0xa00000), // Polyproline II helix, Jmol has no color for it
  T: SecStrColors.betaTurn,
  S: SecStrColors.bend,
});

function makeCarbohydrateSymbolColors(): Record<string, ColorT> {
  const secondaryColor = Color.toHexStyle(SaccharideColors.Secondary);
  const out: Record<string, ColorT> = {};

  for (const [compId, comp] of SaccharideCompIdMap) {
    const primaryColor = Color.toHexStyle(comp.color);
    const hasSecondary = isSaccharideShapeDivided(comp.type);
    out[compId] = (hasSecondary ? `${primaryColor}/${secondaryColor}` : primaryColor) as ColorT;
  }
  return out;
}

export const MvsNamedColorDicts: Record<ColorDictNameT, Record<string, ColorT>> = {
  ElementSymbol: ElementSymbolColorDict,
  ResidueName: ResidueNameColorDict,
  ResidueProperties: ResiduePropertyColorDict,
  SecondaryStructure: SecondaryStructureColors,
  CarbohydrateSymbol: makeCarbohydrateSymbolColors(),
};

export const MvsNamedColorLists: Record<ColorListNameT, ColorT[]> = mapObjectMap(ValuePaletteLists, (colors) =>
  colors.map((color) => Color.toHexStyle(color) as ColorT),
) as Record<ColorListNameT, ColorT[]>;

function MvsColorDict(colors: Record<string, Color>): Record<string, ColorT> {
  return mapObjectMap(colors, (color) => Color.toHexStyle(color) as ColorT);
}
