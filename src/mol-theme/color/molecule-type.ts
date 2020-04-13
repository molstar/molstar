/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorMap } from '../../mol-util/color';
import { StructureElement, Unit, Bond, ElementIndex } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { ColorTheme } from '../color';
import { MoleculeType } from '../../mol-model/structure/model/types';
import { getElementMoleculeType } from '../../mol-model/structure/util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { TableLegend } from '../../mol-util/legend';
import { getAdjustedColorMap } from '../../mol-util/color/color';

export const MoleculeTypeColors = ColorMap({
    water: 0x386cb0,
    ion: 0xf0027f,
    protein: 0xbeaed4,
    RNA: 0xfdc086,
    DNA: 0xbf5b17,
    PNA: 0x42A49A,
    saccharide: 0x7fc97f,
});
export type MoleculeTypeColors = typeof MoleculeTypeColors

const DefaultMoleculeTypeColor = Color(0xffff99);
const Description = 'Assigns a color based on the molecule type of a residue.';

export const MoleculeTypeColorThemeParams = {
    saturation: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 })
};
export type MoleculeTypeColorThemeParams = typeof MoleculeTypeColorThemeParams
export function getMoleculeTypeColorThemeParams(ctx: ThemeDataContext) {
    return MoleculeTypeColorThemeParams; // TODO return copy
}

export function moleculeTypeColor(colorMap: MoleculeTypeColors, unit: Unit, element: ElementIndex): Color {
    const moleculeType = getElementMoleculeType(unit, element);
    switch (moleculeType) {
        case MoleculeType.Water: return colorMap.water;
        case MoleculeType.Ion: return colorMap.ion;
        case MoleculeType.Protein: return colorMap.protein;
        case MoleculeType.RNA: return colorMap.RNA;
        case MoleculeType.DNA: return colorMap.DNA;
        case MoleculeType.PNA: return colorMap.PNA;
        case MoleculeType.Saccharide: return colorMap.saccharide;
    }
    return DefaultMoleculeTypeColor;
}

export function MoleculeTypeColorTheme(ctx: ThemeDataContext, props: PD.Values<MoleculeTypeColorThemeParams>): ColorTheme<MoleculeTypeColorThemeParams> {
    const colorMap = getAdjustedColorMap(MoleculeTypeColors, props.saturation, props.lightness);

    function color(location: Location): Color {
        if (StructureElement.Location.is(location)) {
            return moleculeTypeColor(colorMap, location.unit, location.element);
        } else if (Bond.isLocation(location)) {
            return moleculeTypeColor(colorMap, location.aUnit, location.aUnit.elements[location.aIndex]);
        }
        return DefaultMoleculeTypeColor;
    }

    return {
        factory: MoleculeTypeColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: TableLegend(Object.keys(MoleculeTypeColors).map(name => {
            return [name, (MoleculeTypeColors as any)[name] as Color] as [string, Color];
        }).concat([[ 'Other/unknown', DefaultMoleculeTypeColor ]]))
    };
}

export const MoleculeTypeColorThemeProvider: ColorTheme.Provider<MoleculeTypeColorThemeParams, 'molecule-type'> = {
    name: 'molecule-type',
    label: 'Molecule Type',
    category: ColorTheme.Category.Residue,
    factory: MoleculeTypeColorTheme,
    getParams: getMoleculeTypeColorThemeParams,
    defaultValues: PD.getDefaultValues(MoleculeTypeColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};