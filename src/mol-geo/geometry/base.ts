/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RenderableState } from 'mol-gl/renderable';
import { ValueCell } from 'mol-util';
import { BaseValues } from 'mol-gl/renderable/schema';
import { LocationIterator } from '../util/location-iterator';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { Color } from 'mol-util/color';
import { Vec3 } from 'mol-math/linear-algebra';
import { TransformData, createIdentityTransform } from './transform-data';
import { Theme } from 'mol-theme/theme';
import { ColorNames } from 'mol-util/color/tables';
import { NullLocation } from 'mol-model/location';
import { UniformColorTheme } from 'mol-theme/color/uniform';
import { UniformSizeTheme } from 'mol-theme/size/uniform';

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

//

export namespace BaseGeometry {
    export const Params = {
        alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Opacity' }),
        useFog: PD.Boolean(true),
        highlightColor: PD.Color(Color.fromNormalizedRgb(1.0, 0.4, 0.6)),
        selectColor: PD.Color(Color.fromNormalizedRgb(0.2, 1.0, 0.1)),

        quality: PD.Select<VisualQuality>('auto', VisualQualityOptions),
    }
    export type Params = typeof Params

    export type Counts = { drawCount: number, groupCount: number, instanceCount: number }

    export function createSimple(colorValue = ColorNames.grey, sizeValue = 1, transform?: TransformData) {
        if (!transform) transform = createIdentityTransform()
        const locationIterator = LocationIterator(1, transform.instanceCount.ref.value, () => NullLocation, false, () => false)
        const theme: Theme = {
            color: UniformColorTheme({}, { value: colorValue}),
            size: UniformSizeTheme({}, { value: sizeValue})
        }
        return { transform, locationIterator, theme }
    }

    export function createValues(props: PD.Values<Params>, counts: Counts) {
        return {
            alpha: ValueCell.create(props.alpha),
            uAlpha: ValueCell.create(props.alpha),
            uHighlightColor: ValueCell.create(Color.toArrayNormalized(props.highlightColor, Vec3.zero(), 0)),
            uSelectColor: ValueCell.create(Color.toArrayNormalized(props.selectColor, Vec3.zero(), 0)),
            uGroupCount: ValueCell.create(counts.groupCount),
            drawCount: ValueCell.create(counts.drawCount),
            dUseFog: ValueCell.create(props.useFog),
        }
    }

    export function updateValues(values: BaseValues, props: PD.Values<Params>) {
        if (Color.fromNormalizedArray(values.uHighlightColor.ref.value, 0) !== props.highlightColor) {
            ValueCell.update(values.uHighlightColor, Color.toArrayNormalized(props.highlightColor, values.uHighlightColor.ref.value, 0))
        }
        if (Color.fromNormalizedArray(values.uSelectColor.ref.value, 0) !== props.selectColor) {
            ValueCell.update(values.uSelectColor, Color.toArrayNormalized(props.selectColor, values.uSelectColor.ref.value, 0))
        }
        ValueCell.updateIfChanged(values.alpha, props.alpha) // `uAlpha` is set in renderable.render
        ValueCell.updateIfChanged(values.dUseFog, props.useFog)
    }

    export function createRenderableState(props: Partial<PD.Values<Params>> = {}): RenderableState {
        return {
            visible: true,
            alphaFactor: 1,
            pickable: true,
            opaque: props.alpha === undefined ? true : props.alpha === 1
        }
    }

    export function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
        state.opaque = props.alpha * state.alphaFactor >= 1
    }
}