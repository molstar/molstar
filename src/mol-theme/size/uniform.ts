/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SizeTheme } from '../size';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';

const Description = 'Gives everything the same, uniform size.';

export const UniformSizeThemeParams = {
    value: PD.Numeric(1, { min: 0, max: 20, step: 0.1 }),
};
export type UniformSizeThemeParams = typeof UniformSizeThemeParams
export function getUniformSizeThemeParams(ctx: ThemeDataContext) {
    return UniformSizeThemeParams; // TODO return copy
}

export function UniformSizeTheme(ctx: ThemeDataContext, props: PD.Values<UniformSizeThemeParams>): SizeTheme<UniformSizeThemeParams> {
    const size = props.value;

    return {
        factory: UniformSizeTheme,
        granularity: 'uniform',
        size: () => size,
        props,
        description: Description
    };
}

export const UniformSizeThemeProvider: SizeTheme.Provider<UniformSizeThemeParams, 'uniform'> = {
    name: 'uniform',
    label: 'Uniform',
    category: '',
    factory: UniformSizeTheme,
    getParams: getUniformSizeThemeParams,
    defaultValues: PD.getDefaultValues(UniformSizeThemeParams),
    isApplicable: (ctx: ThemeDataContext) => true
};