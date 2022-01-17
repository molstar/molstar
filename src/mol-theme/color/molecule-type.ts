/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { StructureElement, Unit, Bond, ElementIndex } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { ColorTheme } from '../color';
import { MoleculeType } from '../../mol-model/structure/model/types';
import { getElementMoleculeType } from '../../mol-model/structure/util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { TableLegend } from '../../mol-util/legend';

const DefaultMoleculeTypeColor = Color(0xffff99);
const Description = 'Assigns a color based on the molecule type of a residue.';

export const MoleculeTypeColorThemeParams = {
    saturation: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
    colors: PD.Group({
        water: PD.Color(Color(0x386cb0)),
        ion: PD.Color(Color(0xf0027f)),
        protein: PD.Color(Color(0xbeaed4)),
        rna: PD.Color(Color(0xfdc086)),
        dna: PD.Color(Color(0xbf5b17)),
        pna: PD.Color(Color(0x42A49A)),
        saccharide: PD.Color(Color(0x7fc97f)),
        lipid: PD.Color(Color(0xcccccc)),
    })
};
export type MoleculeTypeColorThemeParams = typeof MoleculeTypeColorThemeParams
export function getMoleculeTypeColorThemeParams(ctx: ThemeDataContext) {
    return MoleculeTypeColorThemeParams; // TODO return copy
}

type MoleculeTypeColorThemeProps = PD.Values<MoleculeTypeColorThemeParams>;
export function moleculeTypeColor(props: MoleculeTypeColorThemeProps, unit: Unit, element: ElementIndex): Color {
    let c = DefaultMoleculeTypeColor;
    const moleculeType = getElementMoleculeType(unit, element);
    switch (moleculeType) {
        case MoleculeType.Water: c = props.colors.water; break;
        case MoleculeType.Ion: c = props.colors.ion; break;
        case MoleculeType.Protein: c = props.colors.protein; break;
        case MoleculeType.RNA: c = props.colors.rna; break;
        case MoleculeType.DNA: c = props.colors.dna; break;
        case MoleculeType.PNA: c = props.colors.pna; break;
        case MoleculeType.Saccharide: c = props.colors.saccharide; break;
        case MoleculeType.Lipid: c = props.colors.lipid; break;
    }
    c = Color.saturate(c, props.saturation);
    c = Color.darken(c, -props.lightness);
    return c;
}

export function MoleculeTypeColorTheme(ctx: ThemeDataContext, props: PD.Values<MoleculeTypeColorThemeParams>): ColorTheme<MoleculeTypeColorThemeParams> {

    function color(location: Location): Color {
        if (StructureElement.Location.is(location)) {
            return moleculeTypeColor(props, location.unit, location.element);
        } else if (Bond.isLocation(location)) {
            return moleculeTypeColor(props, location.aUnit, location.aUnit.elements[location.aIndex]);
        }
        return DefaultMoleculeTypeColor;
    }

    return {
        factory: MoleculeTypeColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: TableLegend(Object.keys(props.colors).map(name => {
            return [name, (props.colors as any)[name] as Color] as [string, Color];
        }).concat([['Other/unknown', DefaultMoleculeTypeColor]]))
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