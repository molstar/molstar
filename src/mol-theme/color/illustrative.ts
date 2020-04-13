/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementSymbol, isNucleic, isProtein, MoleculeType } from '../../mol-model/structure/model/types';
import { Color } from '../../mol-util/color';
import { StructureElement, Unit, Bond } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { ColorTheme } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { elementSymbolColor, ElementSymbolColors } from './element-symbol';
import { getAdjustedColorMap } from '../../mol-util/color/color';

const DefaultIllustrativeColor = Color(0xFFFFFF);
const Description = `Assigns an illustrative color similar to David Goodsell's Molecule of the Month style.`;

export const IllustrativeColorThemeParams = {};
export type IllustrativeColorThemeParams = typeof IllustrativeColorThemeParams
export function getIllustrativeColorThemeParams(ctx: ThemeDataContext) {
    return IllustrativeColorThemeParams; // TODO return copy
}

export function illustrativeColor(colorMap: ElementSymbolColors, typeSymbol: ElementSymbol, moleculeType: MoleculeType) {
    if (isNucleic(moleculeType)) {
        if (typeSymbol === 'O') {
            return Color(0xFF8C8C);
        } else if (typeSymbol === 'P') {
            return Color(0xFF7D7D);
        } else {
            return Color(0xFFA6A6);
        }
    } else if (isProtein(moleculeType)) {
        if (typeSymbol === 'C') {
            return Color(0x7FB2FF);
        } else {
            return Color(0x6699FF);
        }
    } else {
        return elementSymbolColor(colorMap, typeSymbol);
    }
}

export function IllustrativeColorTheme(ctx: ThemeDataContext, props: PD.Values<IllustrativeColorThemeParams>): ColorTheme<IllustrativeColorThemeParams> {
    const colorMap = getAdjustedColorMap(ElementSymbolColors, 0, 0.7);

    function color(location: Location): Color {
        if (StructureElement.Location.is(location)) {
            if (Unit.isAtomic(location.unit)) {
                const moleculeType = location.unit.model.atomicHierarchy.derived.residue.moleculeType[location.unit.residueIndex[location.element]];
                const typeSymbol = location.unit.model.atomicHierarchy.atoms.type_symbol.value(location.element);
                return illustrativeColor(colorMap, typeSymbol, moleculeType);
            }
        } else if (Bond.isLocation(location)) {
            if (Unit.isAtomic(location.aUnit)) {
                const elementIndex = location.aUnit.elements[location.aIndex];
                const moleculeType = location.aUnit.model.atomicHierarchy.derived.residue.moleculeType[location.aUnit.residueIndex[elementIndex]];
                const typeSymbol = location.aUnit.model.atomicHierarchy.atoms.type_symbol.value(elementIndex);
                return illustrativeColor(colorMap, typeSymbol, moleculeType);
            }
        }
        return DefaultIllustrativeColor;
    }

    return {
        factory: IllustrativeColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        // TODO add legend
    };
}

export const IllustrativeColorThemeProvider: ColorTheme.Provider<IllustrativeColorThemeParams, 'illustrative'> = {
    name: 'illustrative',
    label: 'Illustrative',
    category: ColorTheme.Category.Misc,
    factory: IllustrativeColorTheme,
    getParams: getIllustrativeColorThemeParams,
    defaultValues: PD.getDefaultValues(IllustrativeColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};