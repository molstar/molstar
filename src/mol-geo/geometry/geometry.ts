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
import { ColorThemeProps } from 'mol-view/theme/color';
import { LocationIterator } from '../util/location-iterator';
import { ColorType } from './color-data';
import { SizeType } from './size-data';
import { Lines } from './lines/lines';

export type GeometryKindType = { 'mesh': Mesh, 'points': Points, 'lines': Lines }
export type GeometryKind = keyof GeometryKindType
export type Geometry = Helpers.ValueOf<GeometryKindType>

export namespace Geometry {
    export function getDrawCount(geometry: Geometry) {
        switch (geometry.kind) {
            case 'mesh': return geometry.triangleCount * 3
            case 'points': return geometry.pointCount
            case 'lines': return geometry.lineCount
        }
    }

    //

    export const DefaultProps = {
        alpha: 1,
        visible: true,
        depthMask: true,
        useFog: false,
        quality: 'auto' as VisualQuality,
        colorTheme: { name: 'uniform', value: Color(0x22EE11) } as ColorThemeProps,
    }
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

//

export const VisualQualityInfo = {
    'custom': {},
    'auto': {},
    'highest': {},
    'high': {},
    'medium': {},
    'low': {},
    'lowest': {},
}
export type VisualQuality = keyof typeof VisualQualityInfo
export const VisualQualityNames = Object.keys(VisualQualityInfo)