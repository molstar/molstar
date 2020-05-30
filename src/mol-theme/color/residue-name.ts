/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorMap } from '../../mol-util/color';
import { StructureElement, Unit, Bond, ElementIndex } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { ColorTheme } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { TableLegend } from '../../mol-util/legend';
import { getAdjustedColorMap } from '../../mol-util/color/color';

// protein colors from Jmol http://jmol.sourceforge.net/jscolors/
export const ResidueNameColors = ColorMap({
    // standard amino acids
    'ALA': 0x8CFF8C,
    'ARG': 0x00007C,
    'ASN': 0xFF7C70,
    'ASP': 0xA00042,
    'CYS': 0xFFFF70,
    'GLN': 0xFF4C4C,
    'GLU': 0x660000,
    'GLY': 0xEEEEEE,
    'HIS': 0x7070FF,
    'ILE': 0x004C00,
    'LEU': 0x455E45,
    'LYS': 0x4747B8,
    'MET': 0xB8A042,
    'PHE': 0x534C52,
    'PRO': 0x525252,
    'SER': 0xFF7042,
    'THR': 0xB84C00,
    'TRP': 0x4F4600,
    'TYR': 0x8C704C,
    'VAL': 0xFF8CFF,

    // rna bases
    'A': 0xDC143C,  // Crimson Red
    'G': 0x32CD32,  // Lime Green
    'I': 0x9ACD32,  // Yellow Green
    'C': 0xFFD700,  // Gold Yellow
    'T': 0x4169E1,  // Royal Blue
    'U': 0x40E0D0,  // Turquoise Cyan

    // dna bases
    'DA': 0xDC143C,
    'DG': 0x32CD32,
    'DI': 0x9ACD32,
    'DC': 0xFFD700,
    'DT': 0x4169E1,
    'DU': 0x40E0D0,

    // peptide bases
    'APN': 0xDC143C,
    'GPN': 0x32CD32,
    'CPN': 0xFFD700,
    'TPN': 0x4169E1,
});
export type ResidueNameColors = typeof ResidueNameColors

const DefaultResidueNameColor = Color(0xFF00FF);
const Description = 'Assigns a color to every residue according to its name.';

export const ResidueNameColorThemeParams = {
    saturation: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
    lightness: PD.Numeric(1, { min: -6, max: 6, step: 0.1 })
};
export type ResidueNameColorThemeParams = typeof ResidueNameColorThemeParams
export function getResidueNameColorThemeParams(ctx: ThemeDataContext) {
    return ResidueNameColorThemeParams; // TODO return copy
}

function getAtomicCompId(unit: Unit.Atomic, element: ElementIndex) {
    return unit.model.atomicHierarchy.atoms.label_comp_id.value(unit.residueIndex[element]);
}

function getCoarseCompId(unit: Unit.Spheres | Unit.Gaussians, element: ElementIndex) {
    const seqIdBegin = unit.coarseElements.seq_id_begin.value(element);
    const seqIdEnd = unit.coarseElements.seq_id_end.value(element);
    if (seqIdBegin === seqIdEnd) {
        const entityKey = unit.coarseElements.entityKey[element];
        const seq = unit.model.sequence.byEntityKey[entityKey].sequence;
        return seq.compId.value(seqIdBegin - 1); // 1-indexed
    }
}

export function residueNameColor(colorMap: ResidueNameColors, residueName: string): Color {
    const c = colorMap[residueName as keyof ResidueNameColors];
    return c === undefined ? DefaultResidueNameColor : c;
}

export function ResidueNameColorTheme(ctx: ThemeDataContext, props: PD.Values<ResidueNameColorThemeParams>): ColorTheme<ResidueNameColorThemeParams> {
    const colorMap = getAdjustedColorMap(ResidueNameColors, props.saturation, props.lightness);

    function color(location: Location): Color {
        if (StructureElement.Location.is(location)) {
            if (Unit.isAtomic(location.unit)) {
                const compId = getAtomicCompId(location.unit, location.element);
                return residueNameColor(colorMap, compId);
            } else {
                const compId = getCoarseCompId(location.unit, location.element);
                if (compId) return residueNameColor(colorMap, compId);
            }
        } else if (Bond.isLocation(location)) {
            if (Unit.isAtomic(location.aUnit)) {
                const compId = getAtomicCompId(location.aUnit, location.aUnit.elements[location.aIndex]);
                return residueNameColor(colorMap, compId);
            } else {
                const compId = getCoarseCompId(location.aUnit, location.aUnit.elements[location.aIndex]);
                if (compId) return residueNameColor(colorMap, compId);
            }
        }
        return DefaultResidueNameColor;
    }

    return {
        factory: ResidueNameColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: TableLegend(Object.keys(ResidueNameColors).map(name => {
            return [name, (ResidueNameColors as any)[name] as Color] as [string, Color];
        }).concat([[ 'Unknown', DefaultResidueNameColor ]]))
    };
}

export const ResidueNameColorThemeProvider: ColorTheme.Provider<ResidueNameColorThemeParams, 'residue-name'> = {
    name: 'residue-name',
    label: 'Residue Name',
    category: ColorTheme.Category.Residue,
    factory: ResidueNameColorTheme,
    getParams: getResidueNameColorThemeParams,
    defaultValues: PD.getDefaultValues(ResidueNameColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};