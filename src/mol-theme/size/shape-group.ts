/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Location } from '../../mol-model/location';
import { ShapeGroup } from '../../mol-model/shape';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { SizeTheme } from '../../mol-theme/size';

const DefaultSize = 1;
const Description = 'Assigns sizes as defined by the shape object.';

export const ShapeGroupSizeThemeParams = {};
export type ShapeGroupSizeThemeParams = typeof ShapeGroupSizeThemeParams
export function getShapeGroupSizeThemeParams(ctx: ThemeDataContext) {
    return ShapeGroupSizeThemeParams; // TODO return copy
}

export function ShapeGroupSizeTheme(ctx: ThemeDataContext, props: PD.Values<ShapeGroupSizeThemeParams>): SizeTheme<ShapeGroupSizeThemeParams> {
    return {
        factory: ShapeGroupSizeTheme,
        granularity: 'groupInstance',
        size: (location: Location): number => {
            if (ShapeGroup.isLocation(location)) {
                return location.shape.getSize(location.group, location.instance);
            }
            return DefaultSize;
        },
        props,
        description: Description
    };
}

export const ShapeGroupSizeThemeProvider: SizeTheme.Provider<ShapeGroupSizeThemeParams, 'shape-group'> = {
    name: 'shape-group',
    label: 'Shape Group',
    category: '',
    factory: ShapeGroupSizeTheme,
    getParams: getShapeGroupSizeThemeParams,
    defaultValues: PD.getDefaultValues(ShapeGroupSizeThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.shape
};