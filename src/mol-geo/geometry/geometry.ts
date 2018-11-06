/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mesh } from './mesh/mesh';
import { Points } from './points/points';
import { RenderableState } from 'mol-gl/renderable';
import { ValueCell } from 'mol-util';
import { BaseValues } from 'mol-gl/renderable/schema';
import { Color } from 'mol-util/color';
import { ColorThemeOptions, ColorThemeName, ColorScaleOptions, ColorScaleName, ColorTheme } from 'mol-theme/color';
import { LocationIterator } from '../util/location-iterator';
import { ColorType, getColorThemeProps } from './color-data';
import { SizeType, getSizeThemeProps } from './size-data';
import { Lines } from './lines/lines';
import { paramDefaultValues, RangeParam, BooleanParam, SelectParam, ColorParam, NumberParam } from 'mol-util/parameter'
import { DirectVolume } from './direct-volume/direct-volume';
import { SizeTheme, SizeThemeName, SizeThemeOptions } from 'mol-theme/size';

//

export const VisualQualityInfo = {
    'custom': {},
    'auto': {},
    'highest': {},
    'higher': {},
    'high': {},
    'medium': {},
    'low': {},
    'lower': {},
    'lowest': {},
}
export type VisualQuality = keyof typeof VisualQualityInfo
export const VisualQualityNames = Object.keys(VisualQualityInfo)
export const VisualQualityOptions = VisualQualityNames.map(n => [n, n] as [VisualQuality, string])

export interface Theme {
    color: ColorTheme
    size: SizeTheme
}

export function createTheme(props: Geometry.Props) {
    return {
        color: ColorTheme(getColorThemeProps(props)),
        size: SizeTheme(getSizeThemeProps(props))
    }
}

//

export type GeometryKindType = {
    'mesh': Mesh,
    'points': Points,
    'lines': Lines,
    'direct-volume': DirectVolume,
}
export type GeometryKind = keyof GeometryKindType
export type Geometry = Helpers.ValueOf<GeometryKindType>

export namespace Geometry {
    export function getDrawCount(geometry: Geometry) {
        switch (geometry.kind) {
            case 'mesh': return geometry.triangleCount * 3
            case 'points': return geometry.pointCount
            case 'lines': return geometry.lineCount * 2 * 3
            case 'direct-volume': return 12 * 3
        }
    }

    //

    export const Params = {
        alpha: RangeParam('Opacity', '', 1, 0, 1, 0.01),
        visible: BooleanParam('Visible', '', true),
        depthMask: BooleanParam('Depth Mask', '', true),
        useFog: BooleanParam('Use Fog', '', false),
        quality: SelectParam<VisualQuality>('Quality', '', 'auto', VisualQualityOptions),

        colorTheme: SelectParam<ColorThemeName>('Color Theme', '', 'uniform', ColorThemeOptions),
        colorList: SelectParam<ColorScaleName>('Color Scale', '', 'default', ColorScaleOptions),
        colorValue: ColorParam('Color Value', '', Color(0xCCCCCC)),

        sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
        sizeValue: NumberParam('Size Value', '', 1, 0, 20, 0.1),
        sizeFactor: NumberParam('Size Factor', '', 1, 0, 10, 0.1),
    }
    export const DefaultProps = paramDefaultValues(Params)
    export type Props = typeof DefaultProps

    export type Counts = { drawCount: number, groupCount: number, instanceCount: number }

    export function createValues(props: Props, counts: Counts) {
        return {
            uAlpha: ValueCell.create(props.alpha),
            uGroupCount: ValueCell.create(counts.groupCount),
            drawCount: ValueCell.create(counts.drawCount),
            dUseFog: ValueCell.create(props.useFog),
        }
    }

    export function updateValues(values: BaseValues, props: Props) {
        ValueCell.updateIfChanged(values.uAlpha, props.alpha)
        ValueCell.updateIfChanged(values.dUseFog, props.useFog)
    }
}

//

export function createRenderableState(props: Geometry.Props): RenderableState {
    return {
        visible: props.visible,
        depthMask: props.depthMask
    }
}

export function updateRenderableState(state: RenderableState, props: Geometry.Props) {
    state.visible = props.visible
    state.depthMask = props.depthMask
}

//

export function getGranularity(locationIt: LocationIterator, granularity: ColorType | SizeType) {
    // Always use 'group' granularity for 'complex' location iterators,
    // i.e. for which an instance may include multiple units
    return granularity === 'instance' && locationIt.isComplex ? 'group' : granularity
}