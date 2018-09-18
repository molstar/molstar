/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorMap } from 'mol-util/color';
import { StructureElement, Unit, Link, ElementIndex } from 'mol-model/structure';
import { Location } from 'mol-model/location';
import { ColorThemeProps, ColorTheme, TableLegend } from '../color';
import { MoleculeType } from 'mol-model/structure/model/types';
import { getElementMoleculeType } from 'mol-model/structure/util';

const MoleculeTypeColors = ColorMap({
    water: 0x386cb0,
    ion: 0xf0027f,
    protein: 0xbeaed4,
    RNA: 0xfdc086,
    DNA: 0xbf5b17,
    PNA: 0x42A49A,
    saccharide: 0x7fc97f,
})

const DefaultMoleculeTypeColor = Color(0xffff99)
const Description = 'Assigns a color based on the molecule type of a residue.'

export function moleculeTypeColor(unit: Unit, element: ElementIndex): Color {
    const moleculeType = getElementMoleculeType(unit, element)
    switch (moleculeType) {
        case MoleculeType.water: return MoleculeTypeColors.water
        case MoleculeType.ion: return MoleculeTypeColors.ion
        case MoleculeType.protein: return MoleculeTypeColors.protein
        case MoleculeType.RNA: return MoleculeTypeColors.RNA
        case MoleculeType.DNA: return MoleculeTypeColors.DNA
        case MoleculeType.PNA: return MoleculeTypeColors.PNA
        case MoleculeType.saccharide: return MoleculeTypeColors.saccharide
    }
    return DefaultMoleculeTypeColor
}

export function MoleculeTypeColorTheme(props: ColorThemeProps): ColorTheme {
    function color(location: Location): Color {
        if (StructureElement.isLocation(location)) {
            return moleculeTypeColor(location.unit, location.element)
        } else if (Link.isLocation(location)) {
            return moleculeTypeColor(location.aUnit, location.aUnit.elements[location.aIndex])
        }
        return DefaultMoleculeTypeColor
    }

    return {
        granularity: 'group',
        color,
        description: Description,
        legend: TableLegend(Object.keys(MoleculeTypeColors).map(name => {
            return [name, (MoleculeTypeColors as any)[name] as Color] as [string, Color]
        }).concat([[ 'Other/unknown', DefaultMoleculeTypeColor ]]))
    }
}