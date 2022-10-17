/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { StructureElement, Bond, Structure } from '../../mol-model/structure';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { TableLegend, ScaleLegend } from '../../mol-util/legend';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every structure a unique color based on its index.';

export const StructureIndexColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: 'many-distinct' }),
};
export type StructureIndexColorThemeParams = typeof StructureIndexColorThemeParams
export function getStructureIndexColorThemeParams(ctx: ThemeDataContext) {
    return PD.clone(StructureIndexColorThemeParams);
}

export function StructureIndexColorTheme(ctx: ThemeDataContext, props: PD.Values<StructureIndexColorThemeParams>): ColorTheme<StructureIndexColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    if (ctx.structure) {
        const size = (Structure.MaxIndex.get(ctx.structure).value ?? -1) + 1;

        const palette = getPalette(size, props);
        legend = palette.legend;

        color = (location: Location): Color => {
            if (StructureElement.Location.is(location)) {
                return palette.color(Structure.Index.get(location.structure).value || 0)!;
            } else if (Bond.isLocation(location)) {
                return palette.color(Structure.Index.get(location.aStructure).value || 0)!;
            }
            return DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: StructureIndexColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    };
}

export const StructureIndexColorThemeProvider: ColorTheme.Provider<StructureIndexColorThemeParams, 'structure-index'> = {
    name: 'structure-index',
    label: 'Structure Index',
    category: ColorTheme.Category.Chain,
    factory: StructureIndexColorTheme,
    getParams: getStructureIndexColorThemeParams,
    defaultValues: PD.getDefaultValues(StructureIndexColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ctx.structure.elementCount > 0
};