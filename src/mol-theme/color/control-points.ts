/**
 * Copyright (c) 2021-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ColorTheme } from '../color';
import { Color, ColorScale } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { ColorTypeControlPoints } from '../../mol-geo/geometry/color-data';
import { Volume } from '../../mol-model/volume/volume';
import { ColorThemeCategory } from './categories';
import { cpsToColorListRangesEntry as ControlPointsToColorListControlPointsEntry } from '../../mol-plugin-ui/controls/line-graph/line-graph-component';
import { defaultControlPoints } from '../../mol-geo/geometry/direct-volume/direct-volume';

const Description = 'Assign color based on the given value for each volume region defined by control points';

export const ControlPointsThemeName = 'control-points';

export const ControlPointsThemeParams = {
    controlPointsColorList: PD.ColorListControlPoints({
        kind: 'interpolate',
        colors: ControlPointsToColorListControlPointsEntry(defaultControlPoints)
    }, { offsets: true, isEssential: true }),
};
export type ControlPointsThemeParams = typeof ControlPointsThemeParams
export function getControlPointsThemeParams(ctx: ThemeDataContext) {
    return ControlPointsThemeParams; // TODO return copy
}

export function ControlPointsColorTheme(ctx: ThemeDataContext, props: PD.Values<ControlPointsThemeParams>): ColorTheme<ControlPointsThemeParams, ColorTypeControlPoints> {
    const scale = ColorScale.create({ domain: [0, 1], listOrName: props.controlPointsColorList.colors });
    const colors: Color[] = [];
    for (let i = 0; i < 256; ++i) {
        colors[i] = scale.color(i / 255);
    }
    const palette: ColorTheme.Palette = { colors, filter: 'linear' };
    return {
        factory: ControlPointsColorTheme,
        granularity: 'direct',
        props: props,
        description: Description,
        legend: scale.legend,
        palette,
        isRanges: true
    };
}

export const ControlPointsThemeProvider: ColorTheme.Provider<ControlPointsThemeParams, 'control-points'> = {
    name: 'control-points',
    label: 'Control Points',
    category: ColorThemeCategory.Misc,
    factory: ControlPointsColorTheme,
    getParams: getControlPointsThemeParams,
    defaultValues: PD.getDefaultValues(ControlPointsThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.volume && !Volume.Segmentation.get(ctx.volume) && true,
};