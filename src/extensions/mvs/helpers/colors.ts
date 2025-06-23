/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { ElementSymbolColors } from '../../../mol-theme/color/element-symbol';
import { ResidueNameColors } from '../../../mol-theme/color/residue-name';
import { SecondaryStructureColors as SecStrColors } from '../../../mol-theme/color/secondary-structure';
import { Color } from '../../../mol-util/color';
import { ColorList } from '../../../mol-util/color/color';
import { ColorLists } from '../../../mol-util/color/lists';
import { omitObjectKeys } from '../../../mol-util/object';
import { ColorDictNameT, ColorListNameT } from '../tree/mvs/param-types';
import { decodeColor } from './utils';


/** Colors for amino acid groups, based on Clustal (https://www.jalview.org/help/html/colourSchemes/clustal.html) */
const AminoGroupColors = {
    aromatic: decodeColor('#15A4A4')!,
    hydrophobic: decodeColor('#80A0F0')!,
    polar: decodeColor('#15C015')!,
    positive: decodeColor('#F01505')!,
    negative: decodeColor('#C048C0')!,
    proline: decodeColor('#C0C000')!,
    cysteine: decodeColor('#F08080')!,
    glycine: decodeColor('#F09048')!,
};

/** Colors for individual amino acids, based on Clustal (https://www.jalview.org/help/html/colourSchemes/clustal.html), plus Jmol colors for nucleotides (http://jmol.sourceforge.net/jscolors/) */
const ResiduePropertyColors = {
    ...ResidueNameColors,
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
const SecondaryStructureColors = {
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
    P: Color(0xA00000), // Polyproline II helix, Jmol has no color for it
    T: SecStrColors.betaTurn,
    S: SecStrColors.bend,
};

export const MvsNamedColorDicts: Record<ColorDictNameT, Record<string, Color>> = {
    ElementSymbol: omitObjectKeys(ElementSymbolColors, ['C']), // ommitting carbon color to allow easier combination of multiple color layers
    ResidueName: ResidueNameColors,
    ResidueProperties: ResiduePropertyColors,
    SecondaryStructure: SecondaryStructureColors,
};

export const MvsNamedColorLists: Record<ColorListNameT, ColorList> = {
    // Sequential single-hue
    Reds: ColorLists['reds'],
    Oranges: ColorLists['oranges'],
    Greens: ColorLists['greens'],
    Blues: ColorLists['blues'],
    Purples: ColorLists['purples'],
    Greys: ColorLists['greys'],

    // Sequential multi-hue
    OrRd: ColorLists['orange-red'],
    BuGn: ColorLists['blue-green'],
    PuBuGn: ColorLists['purple-blue-green'],
    GnBu: ColorLists['green-blue'],
    PuBu: ColorLists['purple-blue'],
    BuPu: ColorLists['blue-purple'],
    RdPu: ColorLists['red-purple'],
    PuRd: ColorLists['purple-red'],
    YlOrRd: ColorLists['yellow-orange-red'],
    YlOrBr: ColorLists['yellow-orange-brown'],
    YlGn: ColorLists['yellow-green'],
    YlGnBu: ColorLists['yellow-green-blue'],

    Magma: ColorLists['magma'],
    Inferno: ColorLists['inferno'],
    Plasma: ColorLists['plasma'],
    Viridis: ColorLists['viridis'],
    Cividis: ColorLists['cividis'],
    Turbo: ColorLists['turbo'],
    Warm: ColorLists['warm'],
    Cool: ColorLists['cool'],
    CubehelixDefault: ColorLists['cubehelix-default'],

    // Cyclical
    Rainbow: ColorLists['rainbow'],
    Sinebow: ColorLists['sinebow'],

    // Diverging
    RdBu: ColorLists['red-blue'],
    RdGy: ColorLists['red-grey'],
    PiYG: ColorLists['pink-yellow-green'],
    BrBG: ColorLists['brown-white-green'],
    PRGn: ColorLists['purple-green'],
    PuOr: ColorLists['purple-orange'],
    RdYlGn: ColorLists['red-yellow-green'],
    RdYlBu: ColorLists['red-yellow-blue'],
    Spectral: ColorLists['spectral'],

    // Categorical
    Category10: ColorLists['category-10'],
    Observable10: ColorLists['observable-10'],
    Tableau10: ColorLists['tableau-10'],

    Set1: ColorLists['set-1'],
    Set2: ColorLists['set-2'],
    Set3: ColorLists['set-3'],
    Pastel1: ColorLists['pastel-1'],
    Pastel2: ColorLists['pastel-2'],
    Dark2: ColorLists['dark-2'],
    Paired: ColorLists['paired'],
    Accent: ColorLists['accent'],

    // Additional lists, not standard for visualization in general, but commonly used for structures
    Chainbow: ColorLists['turbo-no-black'],
};
