/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme } from './color';
import { SizeTheme } from './size';
import { Structure } from 'mol-model/structure';
import { VolumeData } from 'mol-model/volume';
import { ParamDefinition as PD } from 'mol-util/param-definition';

export interface ThemeRegistryContext {
    colorThemeRegistry: ColorTheme.Registry
    sizeThemeRegistry: SizeTheme.Registry
}

export interface ThemeDataContext {
    [k: string]: any
    structure?: Structure
    volume?: VolumeData
}

export interface Theme {
    color: ColorTheme
    size: SizeTheme
    // label: LabelTheme // TODO
}

type Props = { [k: string]: any }

export function createTheme(ctx: ThemeRegistryContext, data: ThemeDataContext, props: Props, theme?: Theme) {
    theme = theme || createEmptyTheme()

    const colorProps = props.colorTheme as PD.NamedParams
    const sizeProps = props.sizeTheme as PD.NamedParams

    theme.color = ctx.colorThemeRegistry.create(colorProps.name, data, colorProps.params)
    theme.size = ctx.sizeThemeRegistry.create(sizeProps.name, data, sizeProps.params)

    return theme
}

export function createEmptyTheme() {
    return { color: ColorTheme.Empty, size: SizeTheme.Empty }
}