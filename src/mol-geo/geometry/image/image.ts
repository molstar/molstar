/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { hashFnv32a } from '../../../mol-data/util';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { RenderableState } from '../../../mol-gl/renderable';
import { calculateInvariantBoundingSphere, calculateTransformBoundingSphere, TextureImage } from '../../../mol-gl/renderable/util';
import { Sphere3D } from '../../../mol-math/geometry';
import { Vec2, Vec4 } from '../../../mol-math/linear-algebra';
import { Theme } from '../../../mol-theme/theme';
import { ValueCell } from '../../../mol-util';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { BaseGeometry } from '../base';
import { createColors } from '../color-data';
import { GeometryUtils } from '../geometry';
import { createMarkers } from '../marker-data';
import { createEmptyOverpaint } from '../overpaint-data';
import { TransformData } from '../transform-data';
import { createEmptyTransparency } from '../transparency-data';
import { ImageValues } from '../../../mol-gl/renderable/image';
import { fillSerial } from '../../../mol-util/array';
import { createEmptyClipping } from '../clipping-data';

const QuadIndices = new Uint32Array([
    0, 1, 2,
    1, 3, 2
]);

const QuadUvs = new Float32Array([
    0, 1,
    0, 0,
    1, 1,
    1, 0
]);

export const InterpolationTypes = {
    'nearest': 'Nearest',
    'catmulrom': 'Catmulrom (Cubic)',
    'mitchell': 'Mitchell (Cubic)',
    'bspline': 'B-Spline (Cubic)'
};
export type InterpolationTypes = keyof typeof InterpolationTypes;
export const InterpolationTypeNames = Object.keys(InterpolationTypes) as InterpolationTypes[];

export { Image };

interface Image {
    readonly kind: 'image',

    readonly imageTexture: ValueCell<TextureImage<Float32Array>>,
    readonly imageTextureDim: ValueCell<Vec2>,
    readonly cornerBuffer: ValueCell<Float32Array>,
    readonly groupTexture: ValueCell<TextureImage<Float32Array>>,

    /** Bounding sphere of the image */
    boundingSphere: Sphere3D
}

namespace Image {
    export function create(imageTexture: TextureImage<Float32Array>, corners: Float32Array, groupTexture: TextureImage<Float32Array>, image?: Image): Image {
        return image ?
            update(imageTexture, corners, groupTexture, image) :
            fromData(imageTexture, corners, groupTexture);
    }

    function hashCode(image: Image) {
        return hashFnv32a([
            image.cornerBuffer.ref.version
        ]);
    }

    function fromData(imageTexture: TextureImage<Float32Array>, corners: Float32Array, groupTexture: TextureImage<Float32Array>): Image {
        const boundingSphere = Sphere3D();
        let currentHash = -1;

        const width = imageTexture.width;
        const height = imageTexture.height;

        const image = {
            kind: 'image' as const,
            imageTexture: ValueCell.create(imageTexture),
            imageTextureDim: ValueCell.create(Vec2.create(width, height)),
            cornerBuffer: ValueCell.create(corners),
            groupTexture: ValueCell.create(groupTexture),
            get boundingSphere() {
                const newHash = hashCode(image);
                if (newHash !== currentHash) {
                    const b = getBoundingSphere(image.cornerBuffer.ref.value);
                    Sphere3D.copy(boundingSphere, b);
                    currentHash = newHash;
                }
                return boundingSphere;
            },
        };
        return image;
    }

    function update(imageTexture: TextureImage<Float32Array>, corners: Float32Array, groupTexture: TextureImage<Float32Array>, image: Image): Image {
        const width = imageTexture.width;
        const height = imageTexture.height;

        ValueCell.update(image.imageTexture, imageTexture);
        ValueCell.update(image.imageTextureDim, Vec2.set(image.imageTextureDim.ref.value, width, height));
        ValueCell.update(image.cornerBuffer, corners);
        ValueCell.update(image.groupTexture, groupTexture);
        return image;
    }

    export function createEmpty(image?: Image): Image {
        return {} as Image; // TODO
    }

    export const Params = {
        ...BaseGeometry.Params,
        interpolation: PD.Select('bspline', PD.objectToOptions(InterpolationTypes)),
    };
    export type Params = typeof Params

    export const Utils: GeometryUtils<Image, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState,
        updateRenderableState
    };

    function createValues(image: Image, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): ImageValues {

        const { instanceCount, groupCount } = locationIt;
        const color = createColors(locationIt, theme.color);
        const marker = createMarkers(instanceCount * groupCount);
        const overpaint = createEmptyOverpaint();
        const transparency = createEmptyTransparency();
        const clipping = createEmptyClipping();

        const counts = { drawCount: QuadIndices.length, groupCount, instanceCount };

        const invariantBoundingSphere = Sphere3D.clone(image.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform.aTransform.ref.value, instanceCount);

        return {
            ...color,
            ...marker,
            ...overpaint,
            ...transparency,
            ...clipping,
            ...transform,
            ...BaseGeometry.createValues(props, counts),

            aPosition: image.cornerBuffer,
            aUv: ValueCell.create(QuadUvs),
            elements: ValueCell.create(QuadIndices),

            // aGroup is used as a vertex index here, group id is in tGroupTex
            aGroup: ValueCell.create(fillSerial(new Float32Array(4))),
            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),

            dInterpolation: ValueCell.create(props.interpolation),

            uImageTexDim: image.imageTextureDim,
            tImageTex: image.imageTexture,
            tGroupTex: image.groupTexture,
        };
    }

    function createValuesSimple(image: Image, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
        const p = { ...PD.getDefaultValues(Params), ...props };
        return createValues(image, s.transform, s.locationIterator, s.theme, p);
    }

    function updateValues(values: ImageValues, props: PD.Values<Params>) {
        BaseGeometry.updateValues(values, props);
        ValueCell.updateIfChanged(values.dInterpolation, props.interpolation);
    }

    function updateBoundingSphere(values: ImageValues, image: Image) {
        const invariantBoundingSphere = Sphere3D.clone(image.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, values.aTransform.ref.value, values.instanceCount.ref.value);

        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere);
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
            ValueCell.update(values.uInvariantBoundingSphere, Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere));
        }
    }

    function createRenderableState(props: PD.Values<Params>): RenderableState {
        const state = BaseGeometry.createRenderableState(props);
        state.opaque = false;
        return state;
    }

    function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
        BaseGeometry.updateRenderableState(state, props);
        state.opaque = false;
    }
}

//

function getBoundingSphere(corners: Float32Array) {
    return calculateInvariantBoundingSphere(corners, corners.length / 3, 1);
}