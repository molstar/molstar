/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RenderableState } from '../../mol-gl/renderable';
import { ValueCell } from '../../mol-util';
import { BaseValues } from '../../mol-gl/renderable/schema';
import { LocationIterator } from '../util/location-iterator';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { TransformData, createIdentityTransform } from './transform-data';
import { Theme } from '../../mol-theme/theme';
import { ColorNames } from '../../mol-util/color/names';
import { NullLocation } from '../../mol-model/location';
import { UniformColorTheme } from '../../mol-theme/color/uniform';
import { UniformSizeTheme } from '../../mol-theme/size/uniform';

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
};
export type VisualQuality = keyof typeof VisualQualityInfo
export const VisualQualityNames = Object.keys(VisualQualityInfo) as VisualQuality[];
export const VisualQualityOptions = PD.arrayToOptions(VisualQualityNames);

//

export namespace BaseGeometry {
    export const Params = {
        alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Opacity', isEssential: true, description: 'How opaque/transparent the representation is rendered.' }),
        quality: PD.Select<VisualQuality>('auto', VisualQualityOptions, { isEssential: true, description: 'Visual/rendering quality of the representation.' }),
    };
    export type Params = typeof Params

    export const ShadingCategory: PD.Info = { category: 'Shading' };
    export const CustomQualityParamInfo: PD.Info = {
        category: 'Custom Quality',
        hideIf: (params: PD.Values<Params>) => typeof params.quality !== 'undefined' && params.quality !== 'custom'
    };

    export type Counts = { drawCount: number, groupCount: number, instanceCount: number }

    export function createSimple(colorValue = ColorNames.grey, sizeValue = 1, transform?: TransformData) {
        if (!transform) transform = createIdentityTransform();
        const locationIterator = LocationIterator(1, transform.instanceCount.ref.value, () => NullLocation, false, () => false);
        const theme: Theme = {
            color: UniformColorTheme({}, { value: colorValue}),
            size: UniformSizeTheme({}, { value: sizeValue})
        };
        return { transform, locationIterator, theme };
    }

    export function createValues(props: PD.Values<Params>, counts: Counts) {
        return {
            alpha: ValueCell.create(props.alpha),
            uAlpha: ValueCell.create(props.alpha),
            uGroupCount: ValueCell.create(counts.groupCount),
            drawCount: ValueCell.create(counts.drawCount),
        };
    }

    export function updateValues(values: BaseValues, props: PD.Values<Params>) {
        ValueCell.updateIfChanged(values.alpha, props.alpha); // `uAlpha` is set in renderable.render
    }

    export function createRenderableState(props: Partial<PD.Values<Params>> = {}): RenderableState {
        const opaque = props.alpha === undefined ? true : props.alpha === 1;
        return {
            visible: true,
            alphaFactor: 1,
            pickable: true,
            opaque,
            writeDepth: opaque,
        };
    }

    export function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
        state.opaque = props.alpha * state.alphaFactor >= 1;
        state.writeDepth = state.opaque;
    }
}