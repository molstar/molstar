/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SizeTheme } from '../size';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from 'mol-theme/theme';

const Description = 'Gives everything the same, uniform size.'

export const UniformSizeThemeParams = {
    value: PD.Numeric(1, { min: 0, max: 20, step: 0.1 }),
}
export function getUniformSizeThemeParams(ctx: ThemeDataContext) {
    return UniformSizeThemeParams // TODO return copy
}
export type UniformSizeThemeProps = PD.Values<typeof UniformSizeThemeParams>

export function UniformSizeTheme(ctx: ThemeDataContext, props: UniformSizeThemeProps): SizeTheme<UniformSizeThemeProps> {
    const size = props.value

    return {
        granularity: 'uniform',
        size: () => size,
        props,
        description: Description
    }
}

export const UniformSizeThemeProvider: SizeTheme.Provider<typeof UniformSizeThemeParams> = {
    factory: UniformSizeTheme, params: getUniformSizeThemeParams
}