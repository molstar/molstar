/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Location } from '../../mol-model/location';
import { Bond, ElementIndex, Model, StructureElement, Unit } from '../../mol-model/structure';
import { MonosaccharidesColorTable, SaccharideColors, isSaccharideShapeDivided } from '../../mol-model/structure/structure/carbohydrates/constants';
import { Color } from '../../mol-util/color';
import { TableLegend } from '../../mol-util/legend';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import type { ColorTheme, LocationColor } from '../color';
import { ThemeDataContext } from '../theme';
import { ColorThemeCategory } from './categories';


const DefaultColor = Color(0xCCCCCC);
const Description = 'Assigns colors according to the Symbol Nomenclature for Glycans (SNFG).';

export const CarbohydrateSymbolColorThemeParams = {};
export type CarbohydrateSymbolColorThemeParams = typeof CarbohydrateSymbolColorThemeParams
export function getCarbohydrateSymbolColorThemeParams(ctx: ThemeDataContext) {
    return CarbohydrateSymbolColorThemeParams; // TODO return copy
}

export function CarbohydrateSymbolColorTheme(ctx: ThemeDataContext, props: PD.Values<CarbohydrateSymbolColorThemeParams>): ColorTheme<CarbohydrateSymbolColorThemeParams> {
    let color: LocationColor;

    if (ctx.structure) {
        const { elements, getElementIndices } = ctx.structure.carbohydrates;

        const getColor = (unit: Unit, index: ElementIndex, isSecondary: boolean) => {
            if (!Unit.isAtomic(unit)) return DefaultColor;
            const carbs = getElementIndices(unit, index);
            if (carbs.length === 0) return DefaultColor;
            const component = elements[carbs[0]].component;
            if (isSecondary && isSaccharideShapeDivided(component.type)) {
                return SaccharideColors.Secondary;
            } else {
                return component.color;
            }
        };

        color = (location: Location, isSecondary: boolean) => {
            if (StructureElement.Location.is(location)) {
                return getColor(location.unit, location.element, isSecondary);
            } else if (Bond.isLocation(location)) {
                return getColor(location.aUnit, location.aUnit.elements[location.aIndex], isSecondary);
            }
            return DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: CarbohydrateSymbolColorTheme,
        granularity: 'group',
        color: color,
        props: props,
        description: Description,
        legend: TableLegend(MonosaccharidesColorTable)
    };
}

export const CarbohydrateSymbolColorThemeProvider: ColorTheme.Provider<CarbohydrateSymbolColorThemeParams, 'carbohydrate-symbol'> = {
    name: 'carbohydrate-symbol',
    label: 'Carbohydrate Symbol',
    category: ColorThemeCategory.Residue,
    factory: CarbohydrateSymbolColorTheme,
    getParams: getCarbohydrateSymbolColorThemeParams,
    defaultValues: PD.getDefaultValues(CarbohydrateSymbolColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => {
        return !!ctx.structure && ctx.structure.models.some(m => Model.hasCarbohydrate(m));
    }
};