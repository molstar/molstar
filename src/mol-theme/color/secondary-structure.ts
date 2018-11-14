/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorMap } from 'mol-util/color';
import { StructureElement, Unit, Link, ElementIndex } from 'mol-model/structure';
import { Location } from 'mol-model/location';
import { ColorTheme, TableLegend } from '../color';
import { SecondaryStructureType, MoleculeType } from 'mol-model/structure/model/types';
import { getElementMoleculeType } from 'mol-model/structure/util';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from 'mol-theme/theme';

// from Jmol http://jmol.sourceforge.net/jscolors/ (shapely)
const SecondaryStructureColors = ColorMap({
    'alphaHelix': 0xFF0080,
    'threeTenHelix': 0xA00080,
    'piHelix': 0x600080,
    'betaTurn': 0x6080FF,
    'betaStrand': 0xFFC800,
    'coil': 0xFFFFFF,

    'dna': 0xAE00FE,
    'rna': 0xFD0162,

    'carbohydrate': 0xA6A6FA
})

const DefaultSecondaryStructureColor = Color(0x808080)
const Description = 'Assigns a color based on the type of secondary structure and basic molecule type.'

export const SecondaryStructureColorThemeParams = {}
export function getSecondaryStructureColorThemeParams(ctx: ThemeDataContext) {
    return SecondaryStructureColorThemeParams // TODO return copy
}
export type SecondaryStructureColorThemeProps = PD.Values<typeof SecondaryStructureColorThemeParams>

export function secondaryStructureColor(unit: Unit, element: ElementIndex): Color {
    let secStrucType = SecondaryStructureType.create(SecondaryStructureType.Flag.None)
    if (Unit.isAtomic(unit)) {
        secStrucType = unit.model.properties.secondaryStructure.type[unit.residueIndex[element]]
    }

    if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Helix)) {
        if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Helix3Ten)) {
            return SecondaryStructureColors.threeTenHelix
        } else if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.HelixPi)) {
            return SecondaryStructureColors.piHelix
        }
        return SecondaryStructureColors.alphaHelix
    } else if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Beta)) {
        return SecondaryStructureColors.betaStrand
    } else if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Turn)) {
        return SecondaryStructureColors.coil
    } else {
        const moleculeType = getElementMoleculeType(unit, element)
        if (moleculeType === MoleculeType.DNA) {
            return SecondaryStructureColors.dna
        } else if (moleculeType === MoleculeType.RNA) {
            return SecondaryStructureColors.rna
        } else if (moleculeType === MoleculeType.saccharide) {
            return SecondaryStructureColors.carbohydrate
        } else if (moleculeType === MoleculeType.protein) {
            return SecondaryStructureColors.coil
        }
    }
    return DefaultSecondaryStructureColor
}

export function SecondaryStructureColorTheme(ctx: ThemeDataContext, props: SecondaryStructureColorThemeProps): ColorTheme<SecondaryStructureColorThemeProps> {
    function color(location: Location): Color {
        if (StructureElement.isLocation(location)) {
            return secondaryStructureColor(location.unit, location.element)
        } else if (Link.isLocation(location)) {
            return secondaryStructureColor(location.aUnit, location.aUnit.elements[location.aIndex])
        }
        return DefaultSecondaryStructureColor
    }

    return {
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: TableLegend(Object.keys(SecondaryStructureColors).map(name => {
            return [name, (SecondaryStructureColors as any)[name] as Color] as [string, Color]
        }).concat([[ 'Other', DefaultSecondaryStructureColor ]]))
    }
}

export const SecondaryStructureColorThemeProvider: ColorTheme.Provider<typeof SecondaryStructureColorThemeParams> = {
    factory: SecondaryStructureColorTheme, params: getSecondaryStructureColorThemeParams
}