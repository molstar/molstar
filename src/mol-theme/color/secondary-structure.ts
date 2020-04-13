/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorMap } from '../../mol-util/color';
import { StructureElement, Unit, Bond, ElementIndex } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { ColorTheme } from '../color';
import { SecondaryStructureType, MoleculeType } from '../../mol-model/structure/model/types';
import { getElementMoleculeType } from '../../mol-model/structure/util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { TableLegend } from '../../mol-util/legend';
import { SecondaryStructureProvider, SecondaryStructureValue } from '../../mol-model-props/computed/secondary-structure';
import { getAdjustedColorMap } from '../../mol-util/color/color';
import { CustomProperty } from '../../mol-model-props/common/custom-property';

// from Jmol http://jmol.sourceforge.net/jscolors/ (shapely)
const SecondaryStructureColors = ColorMap({
    'alphaHelix': 0xFF0080,
    'threeTenHelix': 0xA00080,
    'piHelix': 0x600080,
    'betaTurn': 0x6080FF,
    'betaStrand': 0xFFC800,
    'coil': 0xFFFFFF,
    'bend': 0x66D8C9 /* biting original color used 0x00FF00 */,
    'turn': 0x00B266,

    'dna': 0xAE00FE,
    'rna': 0xFD0162,

    'carbohydrate': 0xA6A6FA
});
export type SecondaryStructureColors = typeof SecondaryStructureColors

const DefaultSecondaryStructureColor = Color(0x808080);
const Description = 'Assigns a color based on the type of secondary structure and basic molecule type.';

export const SecondaryStructureColorThemeParams = {
    saturation: PD.Numeric(-1, { min: -6, max: 6, step: 0.1 }),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 })
};
export type SecondaryStructureColorThemeParams = typeof SecondaryStructureColorThemeParams
export function getSecondaryStructureColorThemeParams(ctx: ThemeDataContext) {
    return SecondaryStructureColorThemeParams; // TODO return copy
}

export function secondaryStructureColor(colorMap: SecondaryStructureColors, unit: Unit, element: ElementIndex, computedSecondaryStructure?: SecondaryStructureValue): Color {
    let secStrucType = SecondaryStructureType.create(SecondaryStructureType.Flag.None);
    if (computedSecondaryStructure && Unit.isAtomic(unit)) {
        const secondaryStructure = computedSecondaryStructure.get(unit.invariantId);
        if (secondaryStructure) secStrucType = secondaryStructure.type[secondaryStructure.getIndex(unit.residueIndex[element])];
    }

    if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Helix)) {
        if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Helix3Ten)) {
            return colorMap.threeTenHelix;
        } else if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.HelixPi)) {
            return colorMap.piHelix;
        }
        return colorMap.alphaHelix;
    } else if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Beta)) {
        return colorMap.betaStrand;
    } else if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Bend)) {
        return colorMap.bend;
    } else if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Turn)) {
        return colorMap.turn;
    } else {
        const moleculeType = getElementMoleculeType(unit, element);
        if (moleculeType === MoleculeType.DNA) {
            return colorMap.dna;
        } else if (moleculeType === MoleculeType.RNA) {
            return colorMap.rna;
        } else if (moleculeType === MoleculeType.Saccharide) {
            return colorMap.carbohydrate;
        } else if (moleculeType === MoleculeType.Protein) {
            return colorMap.coil;
        }
    }
    return DefaultSecondaryStructureColor;
}

export function SecondaryStructureColorTheme(ctx: ThemeDataContext, props: PD.Values<SecondaryStructureColorThemeParams>): ColorTheme<SecondaryStructureColorThemeParams> {

    const computedSecondaryStructure = ctx.structure && SecondaryStructureProvider.get(ctx.structure);
    const contextHash = computedSecondaryStructure && computedSecondaryStructure.version;

    const colorMap = getAdjustedColorMap(SecondaryStructureColors, props.saturation, props.lightness);

    function color(location: Location): Color {
        if (StructureElement.Location.is(location)) {
            return secondaryStructureColor(colorMap, location.unit, location.element, computedSecondaryStructure?.value);
        } else if (Bond.isLocation(location)) {
            return secondaryStructureColor(colorMap, location.aUnit, location.aUnit.elements[location.aIndex], computedSecondaryStructure?.value);
        }
        return DefaultSecondaryStructureColor;
    }

    return {
        factory: SecondaryStructureColorTheme,
        granularity: 'group',
        color,
        props,
        contextHash,
        description: Description,
        legend: TableLegend(Object.keys(SecondaryStructureColors).map(name => {
            return [name, (SecondaryStructureColors as any)[name] as Color] as [string, Color];
        }).concat([[ 'Other', DefaultSecondaryStructureColor ]]))
    };
}

export const SecondaryStructureColorThemeProvider: ColorTheme.Provider<SecondaryStructureColorThemeParams, 'secondary-structure'> = {
    name: 'secondary-structure',
    label: 'Secondary Structure',
    category: ColorTheme.Category.Residue,
    factory: SecondaryStructureColorTheme,
    getParams: getSecondaryStructureColorThemeParams,
    defaultValues: PD.getDefaultValues(SecondaryStructureColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? SecondaryStructureProvider.attach(ctx, data.structure, void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && data.structure.customPropertyDescriptors.reference(SecondaryStructureProvider.descriptor, false)
    }
};