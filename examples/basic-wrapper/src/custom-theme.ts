import { isPositionLocation } from '@molstar/graphics/geo/util/location-iterator';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { ColorThemeCategory } from '@molstar/graphics/theme/color/categories';
import type { ThemeDataContext } from '@molstar/graphics/theme/theme';
import { Color } from '@molstar/core/util/color';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

export function CustomColorTheme(
    ctx: ThemeDataContext,
    props: PD.Values<{}>
): ColorTheme<{}> {
    const { radius, center } = ctx.structure?.boundary.sphere!;
    const radiusSq = Math.max(radius * radius, 0.001);
    const scale = ColorTheme.PaletteScale;

    return {
        factory: CustomColorTheme,
        granularity: 'vertex',
        color: location => {
            if (!isPositionLocation(location)) return ColorNames.black;
            const dist = Vec3.squaredDistance(location.position, center);
            const t = Math.min(dist / radiusSq, 1);
            return ((t * scale) | 0) as Color;
        },
        palette: {
            filter: 'nearest',
            colors: [
                ColorNames.red,
                ColorNames.pink,
                ColorNames.violet,
                ColorNames.orange,
                ColorNames.yellow,
                ColorNames.green,
                ColorNames.blue
            ]
        },
        props: props,
        description: '',
    };
}

export const CustomColorThemeProvider: ColorTheme.Provider<{}, 'basic-wrapper-custom-color-theme'> = {
    name: 'basic-wrapper-custom-color-theme',
    label: 'Custom Color Theme',
    category: ColorThemeCategory.Misc,
    factory: CustomColorTheme,
    getParams: () => ({}),
    defaultValues: { },
    isApplicable: (ctx: ThemeDataContext) => true,
};
