/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme } from '../color';
import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { isPositionLocation } from '../../mol-geo/util/location-iterator';
import { Vec3 } from '../../mol-math/linear-algebra';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives geometry vertex colors based on positions.';

export const PositionColorThemeParams = {
};
export type PositionColorThemeParams = typeof PositionColorThemeParams
export function getPositionColorThemeParams(ctx: ThemeDataContext) {
    return PositionColorThemeParams; // TODO return copy
}

export function PositionColorTheme(ctx: ThemeDataContext, props: PD.Values<PositionColorThemeParams>): ColorTheme<PositionColorThemeParams> {

    const p = Vec3();
    return {
        factory: PositionColorTheme,
        granularity: 'vertexInstance',
        color: (location: Location): Color => {
            if (isPositionLocation(location)) {
                Vec3.scale(p, location.position, 1 / 50);
                return ColorCosine.palette6(p);
            }
            return DefaultColor;
        },
        props: props,
        description: Description,
    };
}

export const PositionColorThemeProvider: ColorTheme.Provider<PositionColorThemeParams, 'position'> = {
    name: 'position',
    label: 'Position',
    category: ColorTheme.Category.Misc,
    factory: PositionColorTheme,
    getParams: getPositionColorThemeParams,
    defaultValues: PD.getDefaultValues(PositionColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => true
};

//

// ported from https://iquilezles.org/www/articles/palettes/palettes.htm
export namespace ColorCosine {
    const tmp = Vec3();
    export function palette(t: Vec3, a: Vec3, b: Vec3, c: Vec3, d: Vec3): Color {
        // a + b * cos(6.28318 * (c * t + d))
        Vec3.mul(tmp, c, t);
        Vec3.add(tmp, tmp, d);
        Vec3.scale(tmp, tmp, 6.28318);
        tmp[0] = Math.cos(tmp[0]);
        tmp[1] = Math.cos(tmp[1]);
        tmp[2] = Math.cos(tmp[2]);
        Vec3.mul(tmp, b, tmp);
        Vec3.add(tmp, a, tmp);
        return Color.fromNormalizedArray(tmp, 0);
    }

    const a6 = Vec3.create(0.5, 0.5, 0.5);
    const b6 = Vec3.create(0.5, 0.5, 0.5);
    const c6 = Vec3.create(2.0, 1.0, 0.0);
    const d6 = Vec3.create(0.50, 0.20, 0.25);
    export function palette6(t: Vec3) {
        return palette(t, a6, b6, c6, d6);
    }
}