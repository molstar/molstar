/**
 * Copyright (c) 2021-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ColorTheme } from '../color';
import { Color, ColorScale } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { ColorTypeRanges } from '../../mol-geo/geometry/color-data';
import { Volume } from '../../mol-model/volume/volume';
import { ColorThemeCategory } from './categories';
import { generateDefaultControlPoints } from '../../mol-geo/geometry/direct-volume/direct-volume';
import { cpsToColorListRangesEntry } from '../../mol-plugin-ui/controls/line-graph/line-graph-component';

const Description = 'Assign color based on the given value for each volume region defined by control points';

export const defaultControlPoints = generateDefaultControlPoints();

export const RangesColorThemeParams = {
    rangesColorList: PD.ColorListRanges({
        kind: 'interpolate',
        colors: cpsToColorListRangesEntry(defaultControlPoints)
    }, { offsets: true, isEssential: true }),
};
export type RangesColorThemeParams = typeof RangesColorThemeParams
export function getRangesColorThemeParams(ctx: ThemeDataContext) {
    return RangesColorThemeParams; // TODO return copy
}

export function RangesColorTheme(ctx: ThemeDataContext, props: PD.Values<RangesColorThemeParams>): ColorTheme<RangesColorThemeParams, ColorTypeRanges> {
    const scale = ColorScale.create({ domain: [0, 1], listOrName: props.rangesColorList.colors });
    const colors: Color[] = [];
    for (let i = 0; i < 256; ++i) {
        colors[i] = scale.color(i / 255);
    }
    const palette: ColorTheme.Palette = { colors, filter: 'linear' };
    return {
        factory: RangesColorTheme,
        granularity: 'direct',
        props: props,
        description: Description,
        legend: scale.legend,
        palette,
        isRanges: true
    };
}

export const RangesColorThemeProvider: ColorTheme.Provider<RangesColorThemeParams, 'ranges'> = {
    name: 'ranges',
    label: 'Ranges',
    category: ColorThemeCategory.Misc,
    factory: RangesColorTheme,
    getParams: getRangesColorThemeParams,
    defaultValues: PD.getDefaultValues(RangesColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.volume && !Volume.Segmentation.get(ctx.volume) && true,
};