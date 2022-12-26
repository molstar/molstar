/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color, ColorMap } from '../../../../mol-util/color';
import { SaccharideNames } from '../../model/types/saccharides';

// follows community standard from https://www.ncbi.nlm.nih.gov/glycans/snfg.html

export enum SaccharideShape {
    // standard shapes
    FilledSphere, FilledCube, CrossedCube, DividedDiamond, FilledCone, DevidedCone,
    FlatBox, FilledStar, FilledDiamond, FlatDiamond, FlatHexagon, Pentagon,

    // generic shapes for rings with 4, 5, 6, or 7 members
    DiamondPrism, PentagonalPrism, HexagonalPrism, HeptagonalPrism
}

export const SaccharideColors = ColorMap({
    Blue: 0x0090bc,
    Green:	0x00a651,
    Yellow: 0xffd400,
    Orange: 0xf47920,
    Pink: 0xf69ea1,
    Purple: 0xa54399,
    LightBlue: 0x8fcce9,
    Brown: 0xa17a4d,
    Red: 0xed1c24,

    Secondary: 0xf1ece1
});

export enum SaccharideType {
    Hexose, HexNAc, Hexosamine, Hexuronate, Deoxyhexose, DeoxyhexNAc, DiDeoxyhexose,
    Pentose, Deoxynonulosonate, DiDeoxynonulosonate, Unknown, Assigned
}

const SaccharideTypeNameMap = {
    [SaccharideType.Hexose]: 'Hexose',
    [SaccharideType.HexNAc]: 'HexNAc',
    [SaccharideType.Hexosamine]: 'Hexosamine',
    [SaccharideType.Hexuronate]: 'Hexuronate',
    [SaccharideType.Deoxyhexose]: 'Deoxyhexose',
    [SaccharideType.DeoxyhexNAc]: 'DeoxyhexNAc',
    [SaccharideType.DiDeoxyhexose]: 'Di-deoxyhexose',
    [SaccharideType.Pentose]: 'Pentose',
    [SaccharideType.Deoxynonulosonate]: 'Deoxynonulosonate',
    [SaccharideType.DiDeoxynonulosonate]: 'Di-deoxynonulosonate',
    [SaccharideType.Unknown]: 'Unknown',
    [SaccharideType.Assigned]: 'Assigned',
};

export function getSaccharideName(type: SaccharideType) {
    return SaccharideTypeNameMap[type];
}

const SaccharideTypeShapeMap = {
    [SaccharideType.Hexose]: SaccharideShape.FilledSphere,
    [SaccharideType.HexNAc]: SaccharideShape.FilledCube,
    [SaccharideType.Hexosamine]: SaccharideShape.CrossedCube,
    [SaccharideType.Hexuronate]: SaccharideShape.DividedDiamond,
    [SaccharideType.Deoxyhexose]: SaccharideShape.FilledCone,
    [SaccharideType.DeoxyhexNAc]: SaccharideShape.DevidedCone,
    [SaccharideType.DiDeoxyhexose]: SaccharideShape.FlatBox,
    [SaccharideType.Pentose]: SaccharideShape.FilledStar,
    [SaccharideType.Deoxynonulosonate]: SaccharideShape.FilledDiamond,
    [SaccharideType.DiDeoxynonulosonate]: SaccharideShape.FlatDiamond,
    [SaccharideType.Unknown]: SaccharideShape.FlatHexagon,
    [SaccharideType.Assigned]: SaccharideShape.Pentagon,
};

export function getSaccharideShape(type: SaccharideType, ringMemberCount: number): SaccharideShape {
    if (type === SaccharideType.Unknown) {
        if (ringMemberCount === 4) return SaccharideShape.DiamondPrism;
        else if (ringMemberCount === 5) return SaccharideShape.PentagonalPrism;
        else if (ringMemberCount === 6) return SaccharideShape.HexagonalPrism;
        else if (ringMemberCount === 7) return SaccharideShape.HeptagonalPrism;
        else return SaccharideShape.FlatHexagon;
    } else {
        return SaccharideTypeShapeMap[type];
    }
}

export type SaccharideComponent = {
    abbr: string
    name: string
    color: Color
    type: SaccharideType
}

export const UnknownSaccharideComponent: SaccharideComponent = {
    abbr: 'Unk',
    name: 'Unknown',
    color: SaccharideColors.Secondary,
    type: SaccharideType.Unknown
};

const Monosaccharides: SaccharideComponent[] = [
    { abbr: 'Glc', name: 'Glucose', color: SaccharideColors.Blue, type: SaccharideType.Hexose },
    { abbr: 'Man', name: 'Mannose', color: SaccharideColors.Green, type: SaccharideType.Hexose },
    { abbr: 'Gal', name: 'Galactose', color: SaccharideColors.Yellow, type: SaccharideType.Hexose },
    { abbr: 'Gul', name: 'Gulose', color: SaccharideColors.Orange, type: SaccharideType.Hexose },
    { abbr: 'Alt', name: 'Altrose', color: SaccharideColors.Pink, type: SaccharideType.Hexose },
    { abbr: 'All', name: 'Allose', color: SaccharideColors.Purple, type: SaccharideType.Hexose },
    { abbr: 'Tal', name: 'Talose', color: SaccharideColors.LightBlue, type: SaccharideType.Hexose },
    { abbr: 'Ido', name: 'Idose', color: SaccharideColors.Brown, type: SaccharideType.Hexose },

    { abbr: 'GlcNAc', name: 'N-Acetyl Glucosamine', color: SaccharideColors.Blue, type: SaccharideType.HexNAc },
    { abbr: 'ManNAc', name: 'N-Acetyl Mannosamine', color: SaccharideColors.Green, type: SaccharideType.HexNAc },
    { abbr: 'GalNAc', name: 'N-Acetyl Galactosamine', color: SaccharideColors.Yellow, type: SaccharideType.HexNAc },
    { abbr: 'GulNAc', name: 'N-Acetyl Gulosamine', color: SaccharideColors.Orange, type: SaccharideType.HexNAc },
    { abbr: 'AltNAc', name: 'N-Acetyl Altrosamine', color: SaccharideColors.Pink, type: SaccharideType.HexNAc },
    { abbr: 'AllNAc', name: 'N-Acetyl Allosamine', color: SaccharideColors.Purple, type: SaccharideType.HexNAc },
    { abbr: 'TalNAc', name: 'N-Acetyl Talosamine', color: SaccharideColors.LightBlue, type: SaccharideType.HexNAc },
    { abbr: 'IdoNAc', name: 'N-Acetyl Idosamine', color: SaccharideColors.Brown, type: SaccharideType.HexNAc },

    { abbr: 'GlcN', name: 'Glucosamine', color: SaccharideColors.Blue, type: SaccharideType.Hexosamine },
    { abbr: 'ManN', name: 'Mannosamine', color: SaccharideColors.Green, type: SaccharideType.Hexosamine },
    { abbr: 'GalN', name: 'Galactosamine', color: SaccharideColors.Yellow, type: SaccharideType.Hexosamine },
    { abbr: 'GulN', name: 'Gulosamine', color: SaccharideColors.Orange, type: SaccharideType.Hexosamine },
    { abbr: 'AltN', name: 'Altrosamine', color: SaccharideColors.Pink, type: SaccharideType.Hexosamine },
    { abbr: 'AllN', name: 'Allosamine', color: SaccharideColors.Purple, type: SaccharideType.Hexosamine },
    { abbr: 'TalN', name: 'Talosamine', color: SaccharideColors.LightBlue, type: SaccharideType.Hexosamine },
    { abbr: 'IdoN', name: 'Idosamine', color: SaccharideColors.Brown, type: SaccharideType.Hexosamine },

    { abbr: 'GlcA', name: 'Glucuronic Acid', color: SaccharideColors.Blue, type: SaccharideType.Hexuronate },
    { abbr: 'ManA', name: 'Mannuronic Acid', color: SaccharideColors.Green, type: SaccharideType.Hexuronate },
    { abbr: 'GalA', name: 'Galacturonic Acid', color: SaccharideColors.Yellow, type: SaccharideType.Hexuronate },
    { abbr: 'GulA', name: 'Guluronic Acid', color: SaccharideColors.Orange, type: SaccharideType.Hexuronate },
    { abbr: 'AltA', name: 'Altruronic Acid', color: SaccharideColors.Pink, type: SaccharideType.Hexuronate },
    { abbr: 'AllA', name: 'Alluronic Acid', color: SaccharideColors.Purple, type: SaccharideType.Hexuronate },
    { abbr: 'TalA', name: 'Taluronic Acid', color: SaccharideColors.LightBlue, type: SaccharideType.Hexuronate },
    { abbr: 'IdoA', name: 'Iduronic Acid', color: SaccharideColors.Brown, type: SaccharideType.Hexuronate },

    { abbr: 'Qui', name: 'Quinovose', color: SaccharideColors.Blue, type: SaccharideType.Deoxyhexose },
    { abbr: 'Rha', name: 'Rhamnose', color: SaccharideColors.Green, type: SaccharideType.Deoxyhexose },
    { abbr: '6dGul', name: '6-Deoxy Gulose', color: SaccharideColors.Orange, type: SaccharideType.Deoxyhexose },
    { abbr: '6dAlt', name: '6-Deoxy Altrose', color: SaccharideColors.Pink, type: SaccharideType.Deoxyhexose },
    { abbr: '6dTal', name: '6-Deoxy Talose', color: SaccharideColors.LightBlue, type: SaccharideType.Deoxyhexose },
    { abbr: 'Fuc', name: 'Fucose', color: SaccharideColors.Red, type: SaccharideType.Deoxyhexose },

    { abbr: 'QuiNAc', name: 'N-Acetyl Quinovosamine', color: SaccharideColors.Blue, type: SaccharideType.DeoxyhexNAc },
    { abbr: 'RhaNAc', name: 'N-Acetyl Rhamnosamine', color: SaccharideColors.Green, type: SaccharideType.DeoxyhexNAc },
    { abbr: '6dAltNAc', name: 'N-Acetyl 6-Deoxy Altrosamine', color: SaccharideColors.Pink, type: SaccharideType.DeoxyhexNAc },
    { abbr: '6dTalNAc', name: 'N-Acetyl 6-Deoxy Talosamine', color: SaccharideColors.LightBlue, type: SaccharideType.DeoxyhexNAc },
    { abbr: 'FucNAc', name: 'N-Acetyl Fucosamine', color: SaccharideColors.Red, type: SaccharideType.DeoxyhexNAc },

    { abbr: 'Oli', name: 'Olivose', color: SaccharideColors.Blue, type: SaccharideType.DiDeoxyhexose },
    { abbr: 'Tyv', name: 'Tyvelose', color: SaccharideColors.Green, type: SaccharideType.DiDeoxyhexose },
    { abbr: 'Abe', name: 'Abequose', color: SaccharideColors.Orange, type: SaccharideType.DiDeoxyhexose },
    { abbr: 'Par', name: 'Paratose', color: SaccharideColors.Pink, type: SaccharideType.DiDeoxyhexose },
    { abbr: 'Dig', name: 'Digitoxose', color: SaccharideColors.Purple, type: SaccharideType.DiDeoxyhexose },
    { abbr: 'Col', name: 'Colitose', color: SaccharideColors.LightBlue, type: SaccharideType.DiDeoxyhexose },

    { abbr: 'Ara', name: 'Arabinose', color: SaccharideColors.Green, type: SaccharideType.Pentose },
    { abbr: 'Lyx', name: 'Lyxose', color: SaccharideColors.Yellow, type: SaccharideType.Pentose },
    { abbr: 'Xyl', name: 'Xylose', color: SaccharideColors.Orange, type: SaccharideType.Pentose },
    { abbr: 'Rib', name: 'Ribose', color: SaccharideColors.Pink, type: SaccharideType.Pentose },

    { abbr: 'Kdn', name: 'Keto-Deoxy Nonulonic Acid', color: SaccharideColors.Green, type: SaccharideType.Deoxynonulosonate },
    { abbr: 'Neu5Ac', name: 'N-Acetyl Neuraminic Acid', color: SaccharideColors.Purple, type: SaccharideType.Deoxynonulosonate },
    { abbr: 'Neu5Gc', name: 'N-Glycolyl Neuraminic Acid', color: SaccharideColors.LightBlue, type: SaccharideType.Deoxynonulosonate },
    { abbr: 'Neu', name: 'Neuraminic Acid', color: SaccharideColors.Brown, type: SaccharideType.Deoxynonulosonate },
    { abbr: 'Sia', name: 'Sialic acid', color: SaccharideColors.Red, type: SaccharideType.Deoxynonulosonate },

    { abbr: 'Pse', name: 'Pseudaminic Acid', color: SaccharideColors.Green, type: SaccharideType.DiDeoxynonulosonate },
    { abbr: 'Leg', name: 'Legionaminic Acid', color: SaccharideColors.Yellow, type: SaccharideType.DiDeoxynonulosonate },
    { abbr: 'Aci', name: 'Acinetaminic Acid', color: SaccharideColors.Pink, type: SaccharideType.DiDeoxynonulosonate },
    { abbr: '4eLeg', name: '4-Epilegionaminic Acid', color: SaccharideColors.LightBlue, type: SaccharideType.DiDeoxynonulosonate },

    { abbr: 'Bac', name: 'Bacillosamine', color: SaccharideColors.Blue, type: SaccharideType.Unknown },
    { abbr: 'LDmanHep', name: 'L-Glycero-D-Manno Heptose', color: SaccharideColors.Green, type: SaccharideType.Unknown },
    { abbr: 'Kdo', name: 'Keto-Deoxy Octulonic Acid', color: SaccharideColors.Yellow, type: SaccharideType.Unknown },
    { abbr: 'Dha', name: '3-Deoxy Lyxo-Heptulosaric Acid', color: SaccharideColors.Orange, type: SaccharideType.Unknown },
    { abbr: 'DDmanHep', name: 'D-Glycero-D-Manno-Heptose', color: SaccharideColors.Pink, type: SaccharideType.Unknown },
    { abbr: 'MurNAc', name: 'N-Acetyl Muramic Acid', color: SaccharideColors.Purple, type: SaccharideType.Unknown },
    { abbr: 'MurNGc', name: 'N-Glycolyl Muramic Acid', color: SaccharideColors.LightBlue, type: SaccharideType.Unknown },
    { abbr: 'Mur', name: 'Muramic Acid', color: SaccharideColors.Brown, type: SaccharideType.Unknown },

    { abbr: 'Api', name: 'Apicose', color: SaccharideColors.Green, type: SaccharideType.Assigned },
    { abbr: 'Fru', name: 'Fructose', color: SaccharideColors.Green, type: SaccharideType.Assigned },
    { abbr: 'Tag', name: 'Tagatose', color: SaccharideColors.Yellow, type: SaccharideType.Assigned },
    { abbr: 'Sor', name: 'Sorbose', color: SaccharideColors.Orange, type: SaccharideType.Assigned },
    { abbr: 'Psi', name: 'Psicose', color: SaccharideColors.Pink, type: SaccharideType.Assigned },
];

export const SaccharidesSnfgMap = (function () {
    const map = new Map<string, SaccharideComponent>();
    for (let i = 0, il = Monosaccharides.length; i < il; ++i) {
        const saccharide = Monosaccharides[i];
        map.set(saccharide.abbr, saccharide);
    }
    return map;
})();

export const MonosaccharidesColorTable: [string, Color][] = [
    ['Glc-family', SaccharideColors.Blue],
    ['Man-family', SaccharideColors.Green],
    ['Gal-family', SaccharideColors.Yellow],
    ['Gul-family', SaccharideColors.Orange],
    ['Alt-family', SaccharideColors.Pink],
    ['All-family', SaccharideColors.Purple],
    ['Tal-family', SaccharideColors.LightBlue],
    ['Ido-family', SaccharideColors.Brown],
    ['Fuc-family', SaccharideColors.Red],
    ['Generic/Unknown/Secondary', SaccharideColors.Secondary],
];

const CommonSaccharideNames: { [k: string]: string[] } = {
    // Hexose
    Glc: [
        'GLC', 'BGC', 'Z8T',
        'TRE', // di-saccharide but homomer
        'MLR', // tri-saccharide but homomer
    ],
    Man: ['MAN', 'BMA'],
    Gal: ['GLA', 'GAL', 'GZL', 'GXL', 'GIV'],
    Gul: ['4GL', 'GL0', 'GUP', 'Z8H'],
    Alt: ['Z6H', '3MK', 'SHD'],
    All: ['AFD', 'ALL', 'WOO', 'Z2D'],
    Tal: ['ZEE', 'A5C', 'SDY'],
    Ido: ['ZCD', 'Z0F', '4N2'],
    // HexNAc
    GlcNAc: ['NDG', 'NAG', 'NGZ'],
    ManNAc: ['BM3', 'BM7'],
    GalNAc: ['A2G', 'NGA', 'YYQ'],
    GulNAc: ['LXB'],
    AltNAc: [],
    AllNAc: ['NAA'],
    TalNAc: [],
    IdoNAc: ['LXZ', 'HSQ'],
    // Hexosamine
    GlcN: ['PA1', 'GCS'],
    ManN: ['95Z'],
    GalN: ['X6X', '1GN'],
    GulN: [],
    AltN: [],
    AllN: [],
    TalN: [],
    IdoN: [],
    // Hexuronate
    GlcA: ['GCU', 'BDP'],
    ManA: ['MAV', 'BEM'],
    GalA: ['ADA', 'GTR', 'GTK'],
    GulA: ['LGU'],
    AltA: [],
    AllA: [],
    TalA: ['X1X', 'X0X'],
    IdoA: ['IDR'],
    // Deoxyhexose
    Qui: ['G6D', 'YYK'],
    Rha: ['RAM', 'RM4', 'XXR'],
    '6dGul': ['66O'],
    '6dAlt': [],
    '6dTal': [],
    Fuc: ['FUC', 'FUL', 'FCA', 'FCB', 'GYE'],
    // DeoxyhexNAc
    QuiNAc: ['Z9W'],
    RhaNAc: [],
    '6dAltNAc': [],
    '6dTalNAc': [],
    FucNAc: ['49T'],
    // Di-deoxyhexose
    Oli: ['DDA', 'RAE', 'Z5J'],
    Tyv: ['TYV'],
    Abe: ['ABE'],
    Par: ['PZU'],
    Dig: ['Z3U'],
    Col: [],
    // Pentose
    Ara: ['64K', 'ARA', 'ARB', 'AHR', 'FUB', 'BXY', 'BXX', 'SEJ'],
    Lyx: ['LDY', 'Z4W'],
    Xyl: ['XYS', 'XYP', 'XYZ', 'HSY', 'LXC'],
    Rib: ['YYM', 'RIP', 'RIB', 'BDR', '0MK', 'Z6J', '32O'],
    // Deoxynonulosonate
    Kdn: ['KDM', 'KDN'],
    Neu5Ac: ['SIA', 'SLB'],
    Neu5Gc: ['NGC', 'NGE'],
    Neu: [],
    Sia: [],
    // Di-deoxynonulosonate
    Pse: [],
    Leg: [],
    Aci: [],
    '4eLeg': [],
    // Unknown
    Bac: [],
    LDmanHep: ['GMH'],
    Kdo: ['KDO'],
    Dha: [],
    DDmanHep: ['289'],
    MurNAc: ['MUB', 'AMU'],
    MurNGc: [],
    Mur: ['1S4', 'MUR'],
    // Assigned
    Api: ['XXM'],
    Fru: ['BDF', 'Z9N', 'FRU', 'LFR'],
    Tag: ['T6T'],
    Sor: ['SOE', 'UEA'],
    Psi: ['PSV', 'SF6', 'SF9', 'TTV'],
};

/**
 * From http://glycam.org/docs/othertoolsservice/2016/06/09/3d-snfg-list-of-residue-names/#CHARMM
 */
const CharmmSaccharideNames: { [k: string]: string[] } = {
    Glc: ['AGLC', 'BGLC'],
    GlcNAc: ['AGLCNA', 'BGLCNA', 'BGLCN0'],
    GlcA: ['AGLCA', 'BGLCA', 'BGLCA0'],
    Man: ['AMAN', 'BMAN'],
    Rha: ['ARHM', 'BRHM'],
    Ara: ['AARB', 'BARB'],
    Gal: ['AGAL', 'BGAL'],
    GalNAc: ['AGALNA', 'BGALNA'],
    Gul: ['AGUL', 'BGUL'],
    Alt: ['AALT', 'BALT'],
    All: ['AALL', 'BALL'],
    Tal: ['ATAL', 'BTAL'],
    Ido: ['AIDO', 'BIDO'],
    IdoA: ['AIDOA', 'BIDOA'],
    Fuc: ['AFUC', 'BFUC'],
    Lyx: ['ALYF', 'BLYF'],
    Xyl: ['AXYL', 'BXYL', 'AXYF', 'BXYF'],
    Rib: ['ARIB', 'BRIB'],
    Fru: ['AFRU', 'BFRU'],
    Neu5Ac: ['ANE5AC', 'BNE5AC'],
};

/**
 * From http://glycam.org/docs/othertoolsservice/2016/06/09/3d-snfg-list-of-residue-names/#GLYCAM
 */
const GlycamSaccharideNames: { [k: string]: string[] } = {
    Glc: ['0GA', '0GB', '1GA', '1GB', '2GA', '2GB', '3GA', '3GB', '4GA', '4GB', '6GA', '6GB', 'ZGA', 'ZGB', 'YGA', 'YGB', 'XGA', 'XGB', 'WGA', 'WGB', 'VGA', 'VGB', 'UGA', 'UGB', 'TGA', 'TGB', 'SGA', 'SGB', 'RGA', 'RGB', 'QGA', 'QGB', 'PGA', 'PGB', '0gA', '0gB', '1gA', '1gB', '2gA', '2gB', '3gA', '3gB', '4gA', '4gB', '6gA', '6gB', 'ZgA', 'ZgB', 'YgA', 'YgB', 'XgA', 'XgB', 'WgA', 'WgB', 'VgA', 'VgB', 'UgA', 'UgB', 'TgA', 'TgB', 'SgA', 'SgB', 'RgA', 'RgB', 'QgA', 'QgB', 'PgA', 'PgB'],
    GlcNAc: ['0YA', '0YB', '1YA', '1YB', '3YA', '3YB', '4YA', '4YB', '6YA', '6YB', 'WYA', 'WYB', 'VYA', 'VYB', 'UYA', 'UYB', 'QYA', 'QYB', '0yA', '0yB', '1yA', '1yB', '3yA', '3yB', '4yA', '4yB', '6yA', '6yB', 'WyA', 'WyB', 'VyA', 'VyB', 'UyA', 'UyB', 'QyA', 'QyB', '0YS', '0Ys', '3YS', '3Ys', '4YS', '4Ys', '6YS', '6Ys', 'QYS', 'QYs', 'UYS', 'UYs', 'VYS', 'VYs', 'WYS', 'WYs', '0yS', '0ys', '3yS', '3ys', '4yS', '4ys'],
    GlcA: ['0ZA', '0ZB', '1ZA', '1ZB', '2ZA', '2ZB', '3ZA', '3ZB', '4ZA', '4ZB', 'ZZA', 'ZZB', 'YZA', 'YZB', 'WZA', 'WZB', 'TZA', 'TZB', '0zA', '0zB', '1zA', '1zB', '2zA', '2zB', '3zA', '3zB', '4zA', '4zB', 'ZzA', 'ZzB', 'YzA', 'YzB', 'WzA', 'WzB', 'TzA', 'TzB', '0ZBP'],
    GlcN: ['0YN', '3YN', '4YN', '6YN', 'WYN', 'VYN', 'UYN', 'QYN', '3Yn', '4Yn', 'WYn', '0Yn', '0YP', '3YP', '4YP', '6YP', 'WYP', 'VYP', 'UYP', 'QYP', '0Yp', '3Yp', '4Yp', 'WYp'],
    Man: ['0MA', '0MB', '1MA', '1MB', '2MA', '2MB', '3MA', '3MB', '4MA', '4MB', '6MA', '6MB', 'ZMA', 'ZMB', 'YMA', 'YMB', 'XMA', 'XMB', 'WMA', 'WMB', 'VMA', 'VMB', 'UMA', 'UMB', 'TMA', 'TMB', 'SMA', 'SMB', 'RMA', 'RMB', 'QMA', 'QMB', 'PMA', 'PMB', '0mA', '0mB', '1mA', '1mB', '2mA', '2mB', '3mA', '3mB', '4mA', '4mB', '6mA', '6mB', 'ZmA', 'ZmB', 'YmA', 'YmB', 'XmA', 'XmB', 'WmA', 'WmB', 'VmA', 'VmB', 'UmA', 'UmB', 'TmA', 'TmB', 'SmA', 'SmB', 'RmA', 'RmB', 'QmA', 'QmB', 'PmA', 'PmB'],
    ManNAc: ['0WA', '0WB', '1WA', '1WB', '3WA', '3WB', '4WA', '4WB', '6WA', '6WB', 'WWA', 'WWB', 'VWA', 'VWB', 'UWA', 'UWB', 'QWA', 'QWB', '0wA', '0wB', '1wA', '1wB', '3wA', '3wB', '4wA', '4wB', '6wA', '6wB', 'WwA', 'WwB', 'VwA', 'VwB', 'UwA', 'UwB', 'QwA', 'QwB'],
    Ara: ['0AA', '0AB', '1AA', '1AB', '2AA', '2AB', '3AA', '3AB', '4AA', '4AB', 'ZAA', 'ZAB', 'YAA', 'YAB', 'WAA', 'WAB', 'TAA', 'TAB', '0AD', '0AU', '1AD', '1AU', '2AD', '2AU', '3AD', '3AU', '5AD', '5AU', 'ZAD', 'ZAU', '0aA', '0aB', '1aA', '1aB', '2aA', '2aB', '3aA', '3aB', '4aA', '4aB', 'ZaA', 'ZaB', 'YaA', 'YaB', 'WaA', 'WaB', 'TaA', 'TaB', '0aD', '0aU', '1aD', '1aU', '2aD', '2aU', '3aD', '3aU', '5aD', '5aU', 'ZaD', 'ZaU'],
    Gal: ['0LA', '0LB', '1LA', '1LB', '2LA', '2LB', '3LA', '3LB', '4LA', '4LB', '6LA', '6LB', 'ZLA', 'ZLB', 'YLA', 'YLB', 'XLA', 'XLB', 'WLA', 'WLB', 'VLA', 'VLB', 'ULA', 'ULB', 'TLA', 'TLB', 'SLA', 'SLB', 'RLA', 'RLB', 'QLA', 'QLB', 'PLA', 'PLB', '0lA', '0lB', '1lA', '1lB', '2lA', '2lB', '3lA', '3lB', '4lA', '4lB', '6lA', '6lB', 'ZlA', 'ZlB', 'YlA', 'YlB', 'XlA', 'XlB', 'WlA', 'WlB', 'VlA', 'VlB', 'UlA', 'UlB', 'TlA', 'TlB', 'SlA', 'SlB', 'RlA', 'RlB', 'QlA', 'QlB', 'PlA', 'PlB'],
    GalNAc: ['0VA', '0VB', '1VA', '1VB', '3VA', '3VB', '4VA', '4VB', '6VA', '6VB', 'WVA', 'WVB', 'VVA', 'VVB', 'UVA', 'UVB', 'QVA', 'QVB', '0vA', '0vB', '1vA', '1vB', '3vA', '3vB', '4vA', '4vB', '6vA', '6vB', 'WvA', 'WvB', 'VvA', 'VvB', 'UvA', 'UvB', 'QvA', 'QvB'],
    GalA: ['0OA', '0OB', '1OA', '1OB', '2OA', '2OB', '3OA', '3OB', '4OA', '4OB', 'ZOA', 'ZOB', 'YOA', 'YOB', 'WOA', 'WOB', 'TOA', 'TOB', '0oA', '0oB', '1oA', '1oB', '2oA', '2oB', '3oA', '3oB', '4oA', '4oB', 'ZoA', 'ZoB', 'YoA', 'YoB', 'WoA', 'WoB', 'ToA', 'ToB'],
    Gul: ['0KA', '0KB', '1KA', '1KB', '2KA', '2KB', '3KA', '3KB', '4KA', '4KB', '6KA', '6KB', 'ZKA', 'ZKB', 'YKA', 'YKB', 'XKA', 'XKB', 'WKA', 'WKB', 'VKA', 'VKB', 'UKA', 'UKB', 'TKA', 'TKB', 'SKA', 'SKB', 'RKA', 'RKB', 'QKA', 'QKB', 'PKA', 'PKB', '0kA', '0kB', '1kA', '1kB', '2kA', '2kB', '3kA', '3kB', '4kA', '4kB', '6kA', '6kB', 'ZkA', 'ZkB', 'YkA', 'YkB', 'XkA', 'XkB', 'WkA', 'WkB', 'VkA', 'VkB', 'UkA', 'UkB', 'TkA', 'TkB', 'SkA', 'SkB', 'RkA', 'RkB', 'QkA', 'QkB', 'PkA', 'PkB'],
    Alt: ['0EA', '0EB', '1EA', '1EB', '2EA', '2EB', '3EA', '3EB', '4EA', '4EB', '6EA', '6EB', 'ZEA', 'ZEB', 'YEA', 'YEB', 'XEA', 'XEB', 'WEA', 'WEB', 'VEA', 'VEB', 'UEA', 'UEB', 'TEA', 'TEB', 'SEA', 'SEB', 'REA', 'REB', 'QEA', 'QEB', 'PEA', 'PEB', '0eA', '0eB', '1eA', '1eB', '2eA', '2eB', '3eA', '3eB', '4eA', '4eB', '6eA', '6eB', 'ZeA', 'ZeB', 'YeA', 'YeB', 'XeA', 'XeB', 'WeA', 'WeB', 'VeA', 'VeB', 'UeA', 'UeB', 'TeA', 'TeB', 'SeA', 'SeB', 'ReA', 'ReB', 'QeA', 'QeB', 'PeA', 'PeB'],
    All: ['0NA', '0NB', '1NA', '1NB', '2NA', '2NB', '3NA', '3NB', '4NA', '4NB', '6NA', '6NB', 'ZNA', 'ZNB', 'YNA', 'YNB', 'XNA', 'XNB', 'WNA', 'WNB', 'VNA', 'VNB', 'UNA', 'UNB', 'TNA', 'TNB', 'SNA', 'SNB', 'RNA', 'RNB', 'QNA', 'QNB', 'PNA', 'PNB', '0nA', '0nB', '1nA', '1nB', '2nA', '2nB', '3nA', '3nB', '4nA', '4nB', '6nA', '6nB', 'ZnA', 'ZnB', 'YnA', 'YnB', 'XnA', 'XnB', 'WnA', 'WnB', 'VnA', 'VnB', 'UnA', 'UnB', 'TnA', 'TnB', 'SnA', 'SnB', 'RnA', 'RnB', 'QnA', 'QnB', 'PnA', 'PnB'],
    Tal: ['0TA', '0TB', '1TA', '1TB', '2TA', '2TB', '3TA', '3TB', '4TA', '4TB', '6TA', '6TB', 'ZTA', 'ZTB', 'YTA', 'YTB', 'XTA', 'XTB', 'WTA', 'WTB', 'VTA', 'VTB', 'UTA', 'UTB', 'TTA', 'TTB', 'STA', 'STB', 'RTA', 'RTB', 'QTA', 'QTB', 'PTA', 'PTB', '0tA', '0tB', '1tA', '1tB', '2tA', '2tB', '3tA', '3tB', '4tA', '4tB', '6tA', '6tB', 'ZtA', 'ZtB', 'YtA', 'YtB', 'XtA', 'XtB', 'WtA', 'WtB', 'VtA', 'VtB', 'UtA', 'UtB', 'TtA', 'TtB', 'StA', 'StB', 'RtA', 'RtB', 'QtA', 'QtB', 'PtA', 'PtB'],
    Ido: ['0IA', '0IB', '1IA', '1IB', '2IA', '2IB', '3IA', '3IB', '4IA', '4IB', '6IA', '6IB', 'ZIA', 'ZIB', 'YIA', 'YIB', 'XIA', 'XIB', 'WIA', 'WIB', 'VIA', 'VIB', 'UIA', 'UIB', 'TIA', 'TIB', 'SIA', 'SIB', 'RIA', 'RIB', 'QIA', 'QIB', 'PIA', 'PIB', '0iA', '0iB', '1iA', '1iB', '2iA', '2iB', '3iA', '3iB', '4iA', '4iB', '6iA', '6iB', 'ZiA', 'ZiB', 'YiA', 'YiB', 'XiA', 'XiB', 'WiA', 'WiB', 'ViA', 'ViB', 'UiA', 'UiB', 'TiA', 'TiB', 'SiA', 'SiB', 'RiA', 'RiB', 'QiA', 'QiB', 'PiA', 'PiB'],
    IdoA: ['0UA', '0UB', '1UA', '1UB', '2UA', '2UB', '3UA', '3UB', '4UA', '4UB', 'ZUA', 'ZUB', 'YUA', 'YUB', 'WUA', 'WUB', 'TUA', 'TUB', '0uA', '0uB', '1uA', '1uB', '2uA', '2uB', '3uA', '3uB', '4uA', '4uB', 'ZuA', 'ZuB', 'YuA', 'YuB', 'WuA', 'WuB', 'TuA', 'TuB', 'YuAP'],
    Fuc: ['0FA', '0FB', '1FA', '1FB', '2FA', '2FB', '3FA', '3FB', '4FA', '4FB', 'ZFA', 'ZFB', 'YFA', 'YFB', 'WFA', 'WFB', 'TFA', 'TFB', '0fA', '0fB', '1fA', '1fB', '2fA', '2fB', '3fA', '3fB', '4fA', '4fB', 'ZfA', 'ZfB', 'YfA', 'YfB', 'WfA', 'WfB', 'TfA', 'TfB'],
    Rha: ['0HA', '0HB', '1HA', '1HB', '2HA', '2HB', '3HA', '3HB', '4HA', '4HB', 'ZHA', 'ZHB', 'YHA', 'YHB', 'WHA', 'WHB', 'THA', 'THB', '0hA', '0hB', '1hA', '1hB', '2hA', '2hB', '3hA', '3hB', '4hA', '4hB', 'ZhA', 'ZhB', 'YhA', 'YhB', 'WhA', 'WhB', 'ThA', 'ThB'],
    Qui: ['0QA', '0QB', '1QA', '1QB', '2QA', '2QB', '3QA', '3QB', '4QA', '4QB', 'ZQA', 'ZQB', 'YQA', 'YQB', 'WQA', 'WQB', 'TQA', 'TQB', '0qA', '0qB', '1qA', '1qB', '2qA', '2qB', '3qA', '3qB', '4qA', '4qB', 'ZqA', 'ZqB', 'YqA', 'YqB', 'WqA', 'WqB', 'TqA', 'TqB'],
    Lyx: ['0DA', '0DB', '1DA', '1DB', '2DA', '2DB', '3DA', '3DB', '4DA', '4DB', 'ZDA', 'ZDB', 'YDA', 'YDB', 'WDA', 'WDB', 'TDA', 'TDB', '0DD', '0DU', '1DD', '1DU', '2DD', '2DU', '3DD', '3DU', '5DD', '5DU', 'ZDD', 'ZDU', '0dA', '0dB', '1dA', '1dB', '2dA', '2dB', '3dA', '3dB', '4dA', '4dB', 'ZdA', 'ZdB', 'YdA', 'YdB', 'WdA', 'WdB', 'TdA', 'TdB', '0dD', '0dU', '1dD', '1dU', '2dD', '2dU', '3dD', '3dU', '5dD', '5dU', 'ZdD', 'ZdU'],
    Xyl: ['0XA', '0XB', '1XA', '1XB', '2XA', '2XB', '3XA', '3XB', '4XA', '4XB', 'ZXA', 'ZXB', 'YXA', 'YXB', 'WXA', 'WXB', 'TXA', 'TXB', '0XD', '0XU', '1XD', '1XU', '2XD', '2XU', '3XD', '3XU', '5XD', '5XU', 'ZXD', 'ZXU', '0xA', '0xB', '1xA', '1xB', '2xA', '2xB', '3xA', '3xB', '4xA', '4xB', 'ZxA', 'ZxB', 'YxA', 'YxB', 'WxA', 'WxB', 'TxA', 'TxB', '0xD', '0xU', '1xD', '1xU', '2xD', '2xU', '3xD', '3xU', '5xD', '5xU', 'ZxD', 'ZxU'],
    Rib: ['0RA', '0RB', '1RA', '1RB', '2RA', '2RB', '3RA', '3RB', '4RA', '4RB', 'ZRA', 'ZRB', 'YRA', 'YRB', 'WRA', 'WRB', 'TRA', 'TRB', '0RD', '0RU', '1RD', '1RU', '2RD', '2RU', '3RD', '3RU', '5RD', '5RU', 'ZRD', 'ZRU', '0rA', '0rB', '1rA', '1rB', '2rA', '2rB', '3rA', '3rB', '4rA', '4rB', 'ZrA', 'ZrB', 'YrA', 'YrB', 'WrA', 'WrB', 'TrA', 'TrB', '0rD', '0rU', '1rD', '1rU', '2rD', '2rU', '3rD', '3rU', '5rD', '5rU', 'ZrD', 'ZrU'],
    Fru: ['0CA', '0CB', '1CA', '1CB', '2CA', '2CB', '3CA', '3CB', '4CA', '4CB', '5CA', '5CB', 'WCA', 'WCB', '0CD', '0CU', '1CD', '1CU', '2CD', '2CU', '3CD', '3CU', '4CD', '4CU', '6CD', '6CU', 'WCD', 'WCU', 'VCD', 'VCU', 'UCD', 'UCU', 'QCD', 'QCU', '0cA', '0cB', '1cA', '1cB', '2cA', '2cB', '3cA', '3cB', '4cA', '4cB', '5cA', '5cB', 'WcA', 'WcB', '0cD', '0cU', '1cD', '1cU', '2cD', '2cU', '3cD', '3cU', '4cD', '4cU', '6cD', '6cU', 'WcD', 'WcU', 'VcD', 'VcU', 'UcD', 'UcU', 'QcD', 'QcU'],
    Tag: ['0JA', '0JB', '1JA', '1JB', '2JA', '2JB', '3JA', '3JB', '4JA', '4JB', '5JA', '5JB', 'WJA', 'WJB', '0JD', '0JU', '1JD', '1JU', '2JD', '2JU', '3JD', '3JU', '4JD', '4JU', '6JD', '6JU', 'WJD', 'WJU', 'VJD', 'VJU', 'UJD', 'UJU', 'QJD', 'QJU', '0jA', '0jB', '1jA', '1jB', '2jA', '2jB', '3jA', '3jB', '4jA', '4jB', '5jA', '5jB', 'WjA', 'WjB', '0jD', '0jU', '1jD', '1jU', '2jD', '2jU', '3jD', '3jU', '4jD', '4jU', '6jD', '6jU', 'WjD', 'WjU', 'VjD', 'VjU', 'UjD', 'UjU', 'QjD', 'QjU'],
    Sor: ['0BA', '0BB', '1BA', '1BB', '2BA', '2BB', '3BA', '3BB', '4BA', '4BB', '5BA', '5BB', 'WBA', 'WBB', '0BD', '0BU', '1BD', '1BU', '2BD', '2BU', '3BD', '3BU', '4BD', '4BU', '6BD', '6BU', 'WBD', 'WBU', 'VBD', 'VBU', 'UBD', 'UBU', 'QBD', 'QBU', '0bA', '0bB', '1bA', '1bB', '2bA', '2bB', '3bA', '3bB', '4bA', '4bB', '5bA', '5bB', 'WbA', 'WbB', '0bD', '0bU', '1bD', '1bU', '2bD', '2bU', '3bD', '3bU', '4bD', '4bU', '6bD', '6bU', 'WbD', 'WbU', 'VbD', 'VbU', 'UbD', 'UbU', 'QbD', 'QbU'],
    Psi: ['0PA', '0PB', '1PA', '1PB', '2PA', '2PB', '3PA', '3PB', '4PA', '4PB', '5PA', '5PB', 'WPA', 'WPB', '0PD', '0PU', '1PD', '1PU', '2PD', '2PU', '3PD', '3PU', '4PD', '4PU', '6PD', '6PU', 'WPD', 'WPU', 'VPD', 'VPU', 'UPD', 'UPU', 'QPD', 'QPU', '0pA', '0pB', '1pA', '1pB', '2pA', '2pB', '3pA', '3pB', '4pA', '4pB', '5pA', '5pB', 'WpA', 'WpB', '0pD', '0pU', '1pD', '1pU', '2pD', '2pU', '3pD', '3pU', '4pD', '4pU', '6pD', '6pU', 'WpD', 'WpU', 'VpD', 'VpU', 'UpD', 'UpU', 'QpD', 'QpU'],
    Neu5Ac: ['0SA', '0SB', '4SA', '4SB', '7SA', '7SB', '8SA', '8SB', '9SA', '9SB', 'ASA', 'ASB', 'BSA', 'BSB', 'CSA', 'CSB', 'DSA', 'DSB', 'ESA', 'ESB', 'FSA', 'FSB', 'GSA', 'GSB', 'HSA', 'HSB', 'ISA', 'ISB', 'JSA', 'JSB', 'KSA', 'KSB', '0sA', '0sB', '4sA', '4sB', '7sA', '7sB', '8sA', '8sB', '9sA', '9sB', 'AsA', 'AsB', 'BsA', 'BsB', 'CsA', 'CsB', 'DsA', 'DsB', 'EsA', 'EsB', 'FsA', 'FsB', 'GsA', 'GsB', 'HsA', 'HsB', 'IsA', 'IsB', 'JsA', 'JsB', 'KsA', 'KsB'],
    Neu5Gc: ['0GL', '4GL', '7GL', '8GL', '9GL', 'CGL', 'DGL', 'EGL', 'FGL', 'GGL', 'HGL', 'IGL', 'JGL', 'KGL', '0gL', '4gL', '7gL', '8gL', '9gL', 'AgL', 'BgL', 'CgL', 'DgL', 'EgL', 'FgL', 'GgL', 'HgL', 'IgL', 'JgL', 'KgL'],
    Tyv: ['0TV', '0Tv', '1TV', '1Tv', '2TV', '2Tv', '4TV', '4Tv', 'YTV', 'YTv', '0tV', '0tv', '1tV', '1tv', '2tV', '2tv', '4tV', '4tv', 'YtV', 'Ytv'],
    Abe: ['0AE', '2AE', '4AE', 'YGa', '0AF', '2AF', '4AF', 'YAF'],
    Bac: ['0BC', '3BC', '0bC', '3bC'],
    Kdn: ['0KN', '4KN', '5KN', '7KN', '8KN', '9KN', 'AKN', 'BKN', 'CKN', 'DKN', 'EKN', 'FKN', 'GKN', 'HKN', 'IKN', 'JKN', 'KKN', 'LKN', 'MKN', 'NKN', 'OKN', 'PKN', 'QKN', 'RKN', 'SKN', 'TKN', 'UKN', 'VKN', 'WKN', 'XKN', 'YKN', '0Kn', '4Kn', '5Kn', '7Kn', '8Kn', '9Kn', 'AKn', 'BKn', 'CKn', 'DKn', 'EKn', 'FKn', 'GKn', 'HKn', 'IKn', 'JKn', 'KKn', 'LKn', 'MKn', 'NKn', 'OKn', 'PKn', 'QKn', 'RKn', 'SKn', 'TKn', 'UKn', 'VKn', 'WKn', 'XKn', 'YKn'],
    Kdo: ['0KO', '4KO', '5KO', '7KO', '8KO', 'AKO', 'BKO', 'CKO', 'DKO', 'EKO', 'FKO', 'GKO', 'HKO', 'IKO', 'JKO', 'KKO', '0Ko', '4Ko', '5Ko', '7Ko', '8Ko', 'AKo', 'BKo', 'CKo', 'DKo', 'EKo', 'FKo', 'GKo', 'HKo', 'IKo', 'JKo', 'KKo'],
};

const DefaultSaccharideCompIdMap = (function () {
    const map = new Map<string, SaccharideComponent>();
    for (let i = 0, il = Monosaccharides.length; i < il; ++i) {
        const saccharide = Monosaccharides[i];

        const common = CommonSaccharideNames[saccharide.abbr];
        if (common) {
            for (let j = 0, jl = common.length; j < jl; ++j) {
                map.set(common[j], saccharide);
            }
        }

        const charmm = CharmmSaccharideNames[saccharide.abbr];
        if (charmm) {
            for (let j = 0, jl = charmm.length; j < jl; ++j) {
                map.set(charmm[j], saccharide);
            }
        }
    }
    SaccharideNames.forEach(name => {
        if (!map.has(name)) map.set(name, UnknownSaccharideComponent);
    });
    return map;
})();

const GlycamSaccharideCompIdMap = (function () {
    const map = new Map<string, SaccharideComponent>();
    for (let i = 0, il = Monosaccharides.length; i < il; ++i) {
        const saccharide = Monosaccharides[i];

        const common = CommonSaccharideNames[saccharide.abbr];
        if (common) {
            for (let j = 0, jl = common.length; j < jl; ++j) {
                map.set(common[j], saccharide);
            }
        }

        const charmm = CharmmSaccharideNames[saccharide.abbr];
        if (charmm) {
            for (let j = 0, jl = charmm.length; j < jl; ++j) {
                map.set(charmm[j], saccharide);
            }
        }

        const glycam = GlycamSaccharideNames[saccharide.abbr];
        if (glycam) {
            for (let j = 0, jl = glycam.length; j < jl; ++j) {
                // On collision, use PDB name as default.
                if (!map.has(glycam[j])) {
                    map.set(glycam[j], saccharide);
                }
            }
        }
    }
    SaccharideNames.forEach(name => {
        if (!map.has(name)) map.set(name, UnknownSaccharideComponent);
    });
    return map;
})();

export type SaccharideCompIdMapType = 'default' | 'glycam'

export function setSaccharideCompIdMapType(type: SaccharideCompIdMapType) {
    SaccharideCompIdMap = type === 'default' ? DefaultSaccharideCompIdMap : GlycamSaccharideCompIdMap;
}

export let SaccharideCompIdMap = DefaultSaccharideCompIdMap;

export type SaccharideComponentMap = ReadonlyMap<string, SaccharideComponent>
