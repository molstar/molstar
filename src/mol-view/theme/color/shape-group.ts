/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme, ColorThemeProps } from '../color';
import { Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { Shape } from 'mol-model/shape';

const DefaultColor = Color(0xCCCCCC)

export function ShapeGroupColorTheme(props: ColorThemeProps): ColorTheme {
    return {
        granularity: 'group',
        color: (location: Location): Color => {
            if (Shape.isLocation(location)) {
                return location.shape.colors.ref.value[location.group]
            }
            return DefaultColor
        }
    }
}