/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { isSaccharideShapeDivided, SaccharideColors, SaccharideCompIdMap } from '../../../mol-model/structure/structure/carbohydrates/constants';
import { ElementSymbolColors } from '../../../mol-theme/color/element-symbol';
import { ResidueNameColors } from '../../../mol-theme/color/residue-name';
import { SecondaryStructureColors as SecStrColors } from '../../../mol-theme/color/secondary-structure';
import { Color } from '../../../mol-util/color';
import { ColorList } from '../../../mol-util/color/color';
import { ColorLists } from '../../../mol-util/color/lists';
import { mapObjectMap, omitObjectKeys } from '../../../mol-util/object';
import { ColorDictNameT, ColorListNameT, ColorT } from '../tree/mvs/param-types';


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
    P: Color(0xA00000), // Polyproline II helix, Jmol has no color for it
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


export const MvsNamedColorLists: Record<ColorListNameT, ColorT[]> = {
    // Sequential single-hue
    Reds: MvsColorList(ColorLists['reds']),
    Oranges: MvsColorList(ColorLists['oranges']),
    Greens: MvsColorList(ColorLists['greens']),
    Blues: MvsColorList(ColorLists['blues']),
    Purples: MvsColorList(ColorLists['purples']),
    Greys: MvsColorList(ColorLists['greys']),

    // Sequential multi-hue
    OrRd: MvsColorList(ColorLists['orange-red']),
    BuGn: MvsColorList(ColorLists['blue-green']),
    PuBuGn: MvsColorList(ColorLists['purple-blue-green']),
    GnBu: MvsColorList(ColorLists['green-blue']),
    PuBu: MvsColorList(ColorLists['purple-blue']),
    BuPu: MvsColorList(ColorLists['blue-purple']),
    RdPu: MvsColorList(ColorLists['red-purple']),
    PuRd: MvsColorList(ColorLists['purple-red']),
    YlOrRd: MvsColorList(ColorLists['yellow-orange-red']),
    YlOrBr: MvsColorList(ColorLists['yellow-orange-brown']),
    YlGn: MvsColorList(ColorLists['yellow-green']),
    YlGnBu: MvsColorList(ColorLists['yellow-green-blue']),

    Magma: MvsColorList(ColorLists['magma']),
    Inferno: MvsColorList(ColorLists['inferno']),
    Plasma: MvsColorList(ColorLists['plasma']),
    Viridis: MvsColorList(ColorLists['viridis']),
    Cividis: MvsColorList(ColorLists['cividis']),
    Turbo: MvsColorList(ColorLists['turbo']),
    Warm: MvsColorList(ColorLists['warm']),
    Cool: MvsColorList(ColorLists['cool']),
    CubehelixDefault: MvsColorList(ColorLists['cubehelix-default']),

    // Cyclical
    Rainbow: MvsColorList(ColorLists['rainbow']),
    Sinebow: MvsColorList(ColorLists['sinebow']),

    // Diverging
    RdBu: MvsColorList(ColorLists['red-blue']),
    RdGy: MvsColorList(ColorLists['red-grey']),
    PiYG: MvsColorList(ColorLists['pink-yellow-green']),
    BrBG: MvsColorList(ColorLists['brown-white-green']),
    PRGn: MvsColorList(ColorLists['purple-green']),
    PuOr: MvsColorList(ColorLists['purple-orange']),
    RdYlGn: MvsColorList(ColorLists['red-yellow-green']),
    RdYlBu: MvsColorList(ColorLists['red-yellow-blue']),
    Spectral: MvsColorList(ColorLists['spectral']),

    // Categorical
    Category10: MvsColorList(ColorLists['category-10']),
    Observable10: MvsColorList(ColorLists['observable-10']),
    Tableau10: MvsColorList(ColorLists['tableau-10']),

    Set1: MvsColorList(ColorLists['set-1']),
    Set2: MvsColorList(ColorLists['set-2']),
    Set3: MvsColorList(ColorLists['set-3']),
    Pastel1: MvsColorList(ColorLists['pastel-1']),
    Pastel2: MvsColorList(ColorLists['pastel-2']),
    Dark2: MvsColorList(ColorLists['dark-2']),
    Paired: MvsColorList(ColorLists['paired']),
    Accent: MvsColorList(ColorLists['accent']),

    // Additional lists, not standard for visualization in general, but commonly used for structures
    Chainbow: MvsColorList(ColorLists['turbo-no-black']),
};


function MvsColorList(colors: ColorList): ColorT[] {
    return colors.list.map(entry => Color.toHexStyle(Color.fromColorListEntry(entry)) as ColorT);
}

function MvsColorDict(colors: Record<string, Color>): Record<string, ColorT> {
    return mapObjectMap(colors, color => Color.toHexStyle(color) as ColorT);
}
