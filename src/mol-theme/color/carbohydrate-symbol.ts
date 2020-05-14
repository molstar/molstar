/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Bond, ElementIndex, Unit, Model } from '../../mol-model/structure';
import { SaccharideColors, MonosaccharidesColorTable } from '../../mol-model/structure/structure/carbohydrates/constants';
import { Location } from '../../mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { Color } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { TableLegend } from '../../mol-util/legend';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Assigns colors according to the Symbol Nomenclature for Glycans (SNFG).';

export const CarbohydrateSymbolColorThemeParams = { };
export type CarbohydrateSymbolColorThemeParams = typeof CarbohydrateSymbolColorThemeParams
export function getCarbohydrateSymbolColorThemeParams(ctx: ThemeDataContext) {
    return CarbohydrateSymbolColorThemeParams; // TODO return copy
}

export function CarbohydrateSymbolColorTheme(ctx: ThemeDataContext, props: PD.Values<CarbohydrateSymbolColorThemeParams>): ColorTheme<CarbohydrateSymbolColorThemeParams> {
    let color: LocationColor;

    if (ctx.structure) {
        const { elements, getElementIndices } = ctx.structure.carbohydrates;

        const getColor = (unit: Unit, index: ElementIndex) => {
            if (!Unit.isAtomic(unit)) return DefaultColor;
            const carbs = getElementIndices(unit, index);
            return carbs.length > 0 ? elements[carbs[0]].component.color : DefaultColor;
        };

        color = (location: Location, isSecondary: boolean) => {
            if (isSecondary) {
                return SaccharideColors.Secondary;
            } else {
                if (StructureElement.Location.is(location)) {
                    return getColor(location.unit, location.element);
                } else if (Bond.isLocation(location)) {
                    return getColor(location.aUnit, location.aUnit.elements[location.aIndex]);
                }
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
    category: ColorTheme.Category.Residue,
    factory: CarbohydrateSymbolColorTheme,
    getParams: getCarbohydrateSymbolColorThemeParams,
    defaultValues: PD.getDefaultValues(CarbohydrateSymbolColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => {
        return !!ctx.structure && ctx.structure.models.some(m => Model.hasCarbohydrate(m));
    }
};