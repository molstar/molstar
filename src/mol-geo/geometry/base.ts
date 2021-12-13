/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { smoothstep } from '../../mol-math/interpolate';
import { Material } from '../../mol-util/material';

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

export const ColorSmoothingParams = {
    smoothColors: PD.MappedStatic('auto', {
        auto: PD.Group({}),
        on: PD.Group({
            resolutionFactor: PD.Numeric(2, { min: 0.5, max: 6, step: 0.1 }),
            sampleStride: PD.Numeric(3, { min: 1, max: 12, step: 1 }),
        }),
        off: PD.Group({})
    }),
};
export type ColorSmoothingParams = typeof ColorSmoothingParams

export function hasColorSmoothingProp(props: PD.Values<any>): props is PD.Values<ColorSmoothingParams> {
    return !!props.smoothColors;
}

export function getColorSmoothingProps(smoothColors: PD.Values<ColorSmoothingParams>['smoothColors'], preferSmoothing?: boolean, resolution?: number) {
    if ((smoothColors.name === 'on' || (smoothColors.name === 'auto' && preferSmoothing)) && resolution && resolution < 3) {
        let stride = 3;
        if (smoothColors.name === 'on') {
            resolution *= smoothColors.params.resolutionFactor;
            stride = smoothColors.params.sampleStride;
        } else {
            // https://graphtoy.com/?f1(x,t)=(2-smoothstep(0,1.1,x))*x&coords=0.7,0.6,1.8
            resolution *= 2 - smoothstep(0, 1.1, resolution);
            resolution = Math.max(0.5, resolution);
            if (resolution > 1.2) stride = 2;
        }
        return { resolution, stride };
    };
}

//

export namespace BaseGeometry {
    export const MaterialCategory: PD.Info = { category: 'Material' };
    export const ShadingCategory: PD.Info = { category: 'Shading' };
    export const CustomQualityParamInfo: PD.Info = {
        category: 'Custom Quality',
        hideIf: (params: PD.Values<Params>) => typeof params.quality !== 'undefined' && params.quality !== 'custom'
    };

    export const Params = {
        alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Opacity', isEssential: true, description: 'How opaque/transparent the representation is rendered.' }),
        quality: PD.Select<VisualQuality>('auto', VisualQualityOptions, { isEssential: true, description: 'Visual/rendering quality of the representation.' }),
        material: Material.getParam(),
    };
    export type Params = typeof Params

    export type Counts = { drawCount: number, vertexCount: number, groupCount: number, instanceCount: number }

    export function createSimple(colorValue = ColorNames.grey, sizeValue = 1, transform?: TransformData) {
        if (!transform) transform = createIdentityTransform();
        const locationIterator = LocationIterator(1, transform.instanceCount.ref.value, 1, () => NullLocation, false, () => false);
        const theme: Theme = {
            color: UniformColorTheme({}, { value: colorValue }),
            size: UniformSizeTheme({}, { value: sizeValue })
        };
        return { transform, locationIterator, theme };
    }

    export function createValues(props: PD.Values<Params>, counts: Counts) {
        return {
            alpha: ValueCell.create(props.alpha),
            uAlpha: ValueCell.create(props.alpha),
            uVertexCount: ValueCell.create(counts.vertexCount),
            uGroupCount: ValueCell.create(counts.groupCount),
            drawCount: ValueCell.create(counts.drawCount),
            uMetalness: ValueCell.create(props.material.metalness),
            uRoughness: ValueCell.create(props.material.roughness),
            uBumpiness: ValueCell.create(props.material.bumpiness),
            dLightCount: ValueCell.create(1),
        };
    }

    export function updateValues(values: BaseValues, props: PD.Values<Params>) {
        ValueCell.updateIfChanged(values.alpha, props.alpha); // `uAlpha` is set in renderable.render
        ValueCell.updateIfChanged(values.uMetalness, props.material.metalness);
        ValueCell.updateIfChanged(values.uRoughness, props.material.roughness);
        ValueCell.updateIfChanged(values.uBumpiness, props.material.bumpiness);
    }

    export function createRenderableState(props: Partial<PD.Values<Params>> = {}): RenderableState {
        const opaque = props.alpha === undefined ? true : props.alpha === 1;
        return {
            disposed: false,
            visible: true,
            alphaFactor: 1,
            pickable: true,
            colorOnly: false,
            opaque,
            writeDepth: opaque,
            noClip: false,
        };
    }

    export function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
        state.opaque = props.alpha * state.alphaFactor >= 1;
        state.writeDepth = state.opaque;
    }
}