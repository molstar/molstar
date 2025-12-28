/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Lukáš Polák <admin@lukaspolak.cz>
 */

import { Color, ColorMap } from '../../mol-util/color';
import { StructureElement, Unit, Bond, ElementIndex } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import type { ColorTheme } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { TableLegend } from '../../mol-util/legend';
import { getAdjustedColorMap } from '../../mol-util/color/color';
import { getColorMapParams } from '../../mol-util/color/params';
import { ColorThemeCategory } from './categories';

// Colors for charged residues (by-name)
export const ChargedResidueColors = ColorMap({
    // standard amino acids (charged)
    'ARG': 0x0000FF,
    'ASP': 0xFF0000,
    'GLU': 0xFF0000,
    'HIS': 0x33C3F9,
    'LYS': 0x0000FF,

    // standard amino acids (uncharged)
    'ALA': 0xFFFFFF,
    'ASN': 0xFFFFFF,
    'CYS': 0xFFFFFF,
    'GLN': 0xFFFFFF,
    'GLY': 0xFFFFFF,
    'ILE': 0xFFFFFF,
    'LEU': 0xFFFFFF,
    'MET': 0xFFFFFF,
    'PHE': 0xFFFFFF,
    'PRO': 0xFFFFFF,
    'SER': 0xFFFFFF,
    'THR': 0xFFFFFF,
    'TRP': 0xFFFFFF,
    'TYR': 0xFFFFFF,
    'VAL': 0xFFFFFF,

    // common from CCD
    'MSE': 0xFFFFFF,
    'SEP': 0xFFFFFF,
    'TPO': 0xFFFFFF,
    'PTR': 0xFFFFFF,
    'PCA': 0xFFFFFF,
    'HYP': 0xFFFFFF,

    // charmm ff
    'HSD': 0xFFFFFF,
    'HSE': 0xFFFFFF,
    'HSP': 0x0000FF,
    'LSN': 0xFFFFFF,
    'ASPP': 0xFFFFFF,
    'GLUP': 0xFFFFFF,

    // amber ff
    'HID': 0xFFFFFF,
    'HIE': 0xFFFFFF,
    'HIP': 0x0000FF,
    'LYN': 0xFFFFFF,
    'ASH': 0xFFFFFF,
    'GLH': 0xFFFFFF,

    // rna bases
    'A': 0xFFFFFF,
    'G': 0xFFFFFF,
    'I': 0xFFFFFF,
    'C': 0xFFFFFF,
    'T': 0xFFFFFF,
    'U': 0xFFFFFF,

    // dna bases
    'DA': 0xFFFFFF,
    'DG': 0xFFFFFF,
    'DI': 0xFFFFFF,
    'DC': 0xFFFFFF,
    'DT': 0xFFFFFF,
    'DU': 0xFFFFFF,

    // peptide bases
    'APN': 0xFFFFFF,
    'GPN': 0xFFFFFF,
    'CPN': 0xFFFFFF,
    'TPN': 0xFFFFFF,
});
export type ChargedResidueColors = typeof ChargedResidueColors

const DefaultResidueChargeColor = Color(0xFF00FF);
const Description = 'Assigns a color to every residue based on its charge state.';

export const ResidueChargeColorThemeParams = {
    method: PD.MappedStatic('by-name', {
        'by-name': PD.Group({
            saturation: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
            lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
            colors: PD.MappedStatic('default', {
                'default': PD.EmptyGroup(),
                'custom': PD.Group(getColorMapParams(ChargedResidueColors)),
            })
        }, { isFlat: true })
    })
};
export type ResidueChargeColorThemeParams = typeof ResidueChargeColorThemeParams;
export function getResidueChargeColorThemeParams(ctx: ThemeDataContext) {
    return PD.clone(ResidueChargeColorThemeParams);
}

function getAtomicCompId(unit: Unit.Atomic, element: ElementIndex) {
    return unit.model.atomicHierarchy.atoms.label_comp_id.value(element);
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

export function residueChargeColor(colorMap: ColorMap<Record<string, Color>>, residueName: string): Color {
    const c = colorMap[residueName];
    return c === undefined ? DefaultResidueChargeColor : c;
}

export function ResidueChargeColorTheme(ctx: ThemeDataContext, props: PD.Values<ResidueChargeColorThemeParams>): ColorTheme<ResidueChargeColorThemeParams> {
    const { saturation, lightness, colors } = props.method.params;
    const colorMap = getAdjustedColorMap(props.method.params.colors.name === 'default' ? ChargedResidueColors : colors.params, saturation, lightness);

    function color(location: Location): Color {
        if (StructureElement.Location.is(location)) {
            if (Unit.isAtomic(location.unit)) {
                const compId = getAtomicCompId(location.unit, location.element);
                return residueChargeColor(colorMap, compId);
            } else {
                const compId = getCoarseCompId(location.unit, location.element);
                if (compId) return residueChargeColor(colorMap, compId);
            }
        } else if (Bond.isLocation(location)) {
            if (Unit.isAtomic(location.aUnit)) {
                const compId = getAtomicCompId(location.aUnit, location.aUnit.elements[location.aIndex]);
                return residueChargeColor(colorMap, compId);
            } else {
                const compId = getCoarseCompId(location.aUnit, location.aUnit.elements[location.aIndex]);
                if (compId) return residueChargeColor(colorMap, compId);
            }
        }
        return DefaultResidueChargeColor;
    }

    return {
        factory: ResidueChargeColorTheme,
        granularity: 'group',
        preferSmoothing: true,
        color,
        props,
        description: Description,
        legend: TableLegend(Object.keys(colorMap).map(name => {
            return [name, (colorMap as any)[name] as Color] as [string, Color];
        }).concat([['Unknown', DefaultResidueChargeColor]]))
    };
}

export const ResidueChargeColorThemeProvider: ColorTheme.Provider<ResidueChargeColorThemeParams, 'residue-charge'> = {
    name: 'residue-charge',
    label: 'Residue Charge',
    category: ColorThemeCategory.Residue,
    factory: ResidueChargeColorTheme,
    getParams: getResidueChargeColorThemeParams,
    defaultValues: PD.getDefaultValues(ResidueChargeColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};
