/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import type { ColorTheme } from '../../../mol-theme/color';
import { ColorThemeCategory } from '../../../mol-theme/color/categories';
import type { ThemeDataContext } from '../../../mol-theme/theme';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';
import { decodeColor } from '../../../mol-util/color/utils';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { isMVSStructure } from './is-mvs-model-prop';


export function PDSplitColor(info?: PD.Info) {
    return PD.MappedStatic('oneColor', {
        oneColor: PD.Group({
            color: PD.Color(ColorNames.white, { label: ' ' }),
        }, { isFlat: true }),
        splitColor: PD.Group({
            color1: PD.Color(ColorNames.white, { label: '1' }),
            color2: PD.Color(ColorNames.white, { label: '2' }),
        }, { isFlat: true }),
    }, info);
}
export function PDOptionalSplitColor(info?: PD.Info) {
    return PD.MappedStatic('none', {
        none: PD.EmptyGroup(),
        oneColor: PD.Group({
            color: PD.Color(ColorNames.white, { label: ' ' }),
        }, { isFlat: true }),
        splitColor: PD.Group({
            color1: PD.Color(ColorNames.white, { label: '1' }),
            color2: PD.Color(ColorNames.white, { label: '2' }),
        }, { isFlat: true }),
    }, info);
}
export type SplitColorProp = PD.Values<{ x: ReturnType<typeof PDSplitColor> }>['x'];
export type OptionalSplitColorProp = PD.Values<{ x: ReturnType<typeof PDOptionalSplitColor> }>['x'];
type NoneSplitColorProp = Exclude<OptionalSplitColorProp, SplitColorProp>;
type ColorTuple<T extends OptionalSplitColorProp> = T extends SplitColorProp ? [Color, Color] : [undefined, undefined];

const FALLBACK_COLOR = ColorNames.black;
export const SplitColorProp = {
    fromString<T extends string | null>(colorString: T): T extends string ? SplitColorProp : NoneSplitColorProp {
        if (colorString == null) {
            return this.none() satisfies NoneSplitColorProp as any;
        }
        if (colorString.includes('/')) {
            const [c1, c2] = colorString.split('/');
            return { name: 'splitColor', params: { color1: decodeColor(c1) ?? FALLBACK_COLOR, color2: decodeColor(c2) ?? FALLBACK_COLOR } } satisfies SplitColorProp as any;
        } else {
            return { name: 'oneColor', params: { color: decodeColor(colorString) ?? FALLBACK_COLOR } } satisfies SplitColorProp as any;
        }
    },
    none(): NoneSplitColorProp {
        return { name: 'none', params: {} };
    },
    toString(value: OptionalSplitColorProp): string {
        if (value.name === 'none') return 'None';
        if (value.name === 'oneColor') return Color.toHexStyle(value.params.color);
        else return `${Color.toHexStyle(value.params.color1)}/${Color.toHexStyle(value.params.color2)}`;
    },
    toTuple<T extends OptionalSplitColorProp>(value: T): ColorTuple<T> {
        if (value.name === 'none') return [undefined, undefined] as ColorTuple<T>;
        if (value.name === 'oneColor') return [value.params.color, value.params.color] as ColorTuple<T>;
        else return [value.params.color1, value.params.color2] as ColorTuple<T>;
    },
};


export const MVSSplitUniformColorThemeParams = {
    value: PDSplitColor(),
};
export type MVSSplitUniformColorThemeParams = typeof MVSSplitUniformColorThemeParams;
export type MVSSplitUniformColorThemeProps = PD.Values<MVSSplitUniformColorThemeParams>;

export function MVSSplitUniformColorTheme(ctx: ThemeDataContext, props: PD.Values<MVSSplitUniformColorThemeParams>): ColorTheme<MVSSplitUniformColorThemeParams> {
    const [color1, color2] = SplitColorProp.toTuple(props.value);

    return {
        factory: MVSSplitUniformColorTheme,
        granularity: 'groupInstance',
        color: (location, isSecondary) => isSecondary ? color2 : color1,
        props: props,
        description: 'Gives everything the same, uniform color, optionally split into two colors (when representation supports split coloring, e.g. for carbohydrate symbols, nucleic acid cartoon)',
    };
}

export const MVSSplitUniformColorThemeProvider: ColorTheme.Provider<MVSSplitUniformColorThemeParams, 'mvs-split-uniform'> = {
    name: 'mvs-split-uniform',
    label: 'MVS Split Uniform',
    category: ColorThemeCategory.Misc,
    factory: MVSSplitUniformColorTheme,
    getParams: () => MVSSplitUniformColorThemeParams,
    defaultValues: PD.getDefaultValues(MVSSplitUniformColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && isMVSStructure(ctx.structure),
};

