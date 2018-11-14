/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme } from './color';
import { SizeTheme } from './size';
import { Structure } from 'mol-model/structure';
import { VolumeData } from 'mol-model/volume';

export interface ThemeRegistryContext {
    colorThemeRegistry: ColorTheme.Registry
    sizeThemeRegistry: SizeTheme.Registry
}

export interface ThemeDataContext {
    [k: string]: any
    structure?: Structure
    volume?: VolumeData
}

export interface ThemeProps { color?: {}, size?: {} }

export interface Theme {
    color: ColorTheme
    size: SizeTheme
    // label: LabelTheme // TODO
}

type Props = { [k: string]: any }

export function createTheme(ctx: ThemeRegistryContext, data: ThemeDataContext, props: Props, themeProps: ThemeProps = {}, theme?: Theme) {
    theme = theme || {
        color: ColorTheme.Empty,
        size: SizeTheme.Empty
    }
    // TODO check if props have changed
    if (typeof props.colorTheme === 'string') {
        theme.color = ctx.colorThemeRegistry.create(props.colorTheme, data, themeProps.color)
    }
    if (typeof props.sizeTheme === 'string') {
        theme.size = ctx.sizeThemeRegistry.create(props.sizeTheme, data, themeProps.size)
    }
    return theme
}