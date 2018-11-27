/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme } from '../color';
import { Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { Shape } from 'mol-model/shape';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from 'mol-theme/theme';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Assigns colors as defined by the shape object.'

export const ShapeGroupColorThemeParams = {}
export type ShapeGroupColorThemeParams = typeof ShapeGroupColorThemeParams
export function getShapeGroupColorThemeParams(ctx: ThemeDataContext) {
    return ShapeGroupColorThemeParams // TODO return copy
}

export function ShapeGroupColorTheme(ctx: ThemeDataContext, props: PD.Values<ShapeGroupColorThemeParams>): ColorTheme<ShapeGroupColorThemeParams> {
    return {
        factory: ShapeGroupColorTheme,
        granularity: 'group',
        color: (location: Location): Color => {
            if (Shape.isLocation(location)) {
                return location.shape.colors.ref.value[location.group]
            }
            return DefaultColor
        },
        props,
        description: Description
    }
}

export const ShapeGroupColorThemeProvider: ColorTheme.Provider<ShapeGroupColorThemeParams> = {
    label: 'Shape Group',
    factory: ShapeGroupColorTheme,
    getParams: getShapeGroupColorThemeParams
}