/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { StructureElement, Bond } from '../../mol-model/structure';
import { OrderedSet } from '../../mol-data/int';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { TableLegend, ScaleLegend } from '../../mol-util/legend';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every element (atom or coarse sphere/gaussian) a unique color based on the position (index) of the element in the list of elements in the structure.';

export const ElementIndexColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: 'red-yellow-blue' }),
};
export type ElementIndexColorThemeParams = typeof ElementIndexColorThemeParams
export function getElementIndexColorThemeParams(ctx: ThemeDataContext) {
    return ElementIndexColorThemeParams; // TODO return copy
}

export function ElementIndexColorTheme(ctx: ThemeDataContext, props: PD.Values<ElementIndexColorThemeParams>): ColorTheme<ElementIndexColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    if (ctx.structure) {
        const { units } = ctx.structure.root;
        const unitCount = units.length;
        const cummulativeElementCount = new Map<number, number>();
        const unitIdIndex = new Map<number, number>();

        let elementCount = 0;
        for (let i = 0; i < unitCount; ++i) {
            cummulativeElementCount.set(i, elementCount);
            elementCount += units[i].elements.length;
            unitIdIndex.set(units[i].id, i);
        }

        const palette = getPalette(elementCount, props);
        legend = palette.legend;

        color = (location: Location): Color => {
            if (StructureElement.Location.is(location)) {
                const unitIndex = unitIdIndex.get(location.unit.id)!;
                const unitElementIndex = OrderedSet.findPredecessorIndex(units[unitIndex].elements, location.element);
                return palette.color(cummulativeElementCount.get(unitIndex)! + unitElementIndex);
            } else if (Bond.isLocation(location)) {
                const unitIndex = unitIdIndex.get(location.aUnit.id)!;
                const unitElementIndex = OrderedSet.findPredecessorIndex(units[unitIndex].elements, location.aUnit.elements[location.aIndex]);
                return palette.color(cummulativeElementCount.get(unitIndex)! + unitElementIndex);
            }
            return DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: ElementIndexColorTheme,
        granularity: 'groupInstance',
        color,
        props,
        description: Description,
        legend
    };
}

export const ElementIndexColorThemeProvider: ColorTheme.Provider<ElementIndexColorThemeParams, 'element-index'> = {
    name: 'element-index',
    label: 'Element Index',
    category: ColorTheme.Category.Atom,
    factory: ElementIndexColorTheme,
    getParams: getElementIndexColorThemeParams,
    defaultValues: PD.getDefaultValues(ElementIndexColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};