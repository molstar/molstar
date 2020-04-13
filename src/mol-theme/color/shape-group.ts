/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme } from '../color';
import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { ShapeGroup } from '../../mol-model/shape';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Assigns colors as defined by the shape object.';

export const ShapeGroupColorThemeParams = {};
export type ShapeGroupColorThemeParams = typeof ShapeGroupColorThemeParams
export function getShapeGroupColorThemeParams(ctx: ThemeDataContext) {
    return ShapeGroupColorThemeParams; // TODO return copy
}

export function ShapeGroupColorTheme(ctx: ThemeDataContext, props: PD.Values<ShapeGroupColorThemeParams>): ColorTheme<ShapeGroupColorThemeParams> {
    return {
        factory: ShapeGroupColorTheme,
        granularity: 'groupInstance',
        color: (location: Location): Color => {
            if (ShapeGroup.isLocation(location)) {
                return location.shape.getColor(location.group, location.instance);
            }
            return DefaultColor;
        },
        props,
        description: Description
    };
}

export const ShapeGroupColorThemeProvider: ColorTheme.Provider<ShapeGroupColorThemeParams, 'shape-group'> = {
    name: 'shape-group',
    label: 'Shape Group',
    category: ColorTheme.Category.Misc,
    factory: ShapeGroupColorTheme,
    getParams: getShapeGroupColorThemeParams,
    defaultValues: PD.getDefaultValues(ShapeGroupColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.shape
};