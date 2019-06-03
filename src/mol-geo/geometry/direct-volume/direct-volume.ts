/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../../mol-util'
import { Sphere3D, Box3D } from '../../../mol-math/geometry'
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { DirectVolumeValues } from '../../../mol-gl/renderable/direct-volume';
import { Vec3, Mat4, Vec2 } from '../../../mol-math/linear-algebra';
import { Box } from '../../primitive/box';
import { createTransferFunctionTexture, getControlPointsFromVec2Array } from './transfer-function';
import { Texture } from '../../../mol-gl/webgl/texture';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { TransformData } from '../transform-data';
import { createColors } from '../color-data';
import { createMarkers } from '../marker-data';
import { GeometryUtils } from '../geometry';
import { transformPositionArray } from '../../../mol-geo/util';
import { calculateBoundingSphere } from '../../../mol-gl/renderable/util';
import { Theme } from '../../../mol-theme/theme';
import { RenderableState } from '../../../mol-gl/renderable';
import { ColorListOptions, ColorListName } from '../../../mol-util/color/scale';
import { Color } from '../../../mol-util/color';
import { BaseGeometry } from '../base';
import { createEmptyOverpaint } from '../overpaint-data';
import { createEmptyTransparency } from '../transparency-data';

const VolumeBox = Box()
const RenderModeOptions = [['isosurface', 'Isosurface'], ['volume', 'Volume']] as [string, string][]

export interface DirectVolume {
    readonly kind: 'direct-volume',
    readonly gridTexture: ValueCell<Texture>,
    readonly gridTextureDim: ValueCell<Vec3>,
    readonly gridDimension: ValueCell<Vec3>,
    readonly bboxSize: ValueCell<Vec3>
    readonly bboxMin: ValueCell<Vec3>
    readonly bboxMax: ValueCell<Vec3>
    readonly transform: ValueCell<Mat4>

    /** Bounding sphere of the volume */
    boundingSphere?: Sphere3D
}

export namespace DirectVolume {
    export function create(bbox: Box3D, gridDimension: Vec3, transform: Mat4, texture: Texture, directVolume?: DirectVolume): DirectVolume {
        const { width, height, depth } = texture
        if (directVolume) {
            ValueCell.update(directVolume.gridDimension, gridDimension)
            ValueCell.update(directVolume.gridTextureDim, Vec3.set(directVolume.gridTextureDim.ref.value, width, height, depth))
            ValueCell.update(directVolume.bboxMin, bbox.min)
            ValueCell.update(directVolume.bboxMax, bbox.max)
            ValueCell.update(directVolume.bboxSize, Vec3.sub(directVolume.bboxSize.ref.value, bbox.max, bbox.min))
            ValueCell.update(directVolume.transform, transform)
            return directVolume
        } else {
            return {
                kind: 'direct-volume',
                gridDimension: ValueCell.create(gridDimension),
                gridTexture: ValueCell.create(texture),
                gridTextureDim: ValueCell.create(Vec3.create(width, height, depth)),
                bboxMin: ValueCell.create(bbox.min),
                bboxMax: ValueCell.create(bbox.max),
                bboxSize: ValueCell.create(Vec3.sub(Vec3.zero(), bbox.max, bbox.min)),
                transform: ValueCell.create(transform),
            }
        }
    }

    export function createEmpty(directVolume?: DirectVolume): DirectVolume {
        return {} as DirectVolume // TODO
    }

    export const Params = {
        ...BaseGeometry.Params,
        isoValueNorm: PD.Numeric(0.22, { min: 0, max: 1, step: 0.01 }, { description: 'Normalized Isolevel Value' }),
        renderMode: PD.Select('volume', RenderModeOptions),
        controlPoints: PD.LineGraph([
            Vec2.create(0.19, 0.0), Vec2.create(0.2, 0.05), Vec2.create(0.25, 0.05), Vec2.create(0.26, 0.0),
            Vec2.create(0.79, 0.0), Vec2.create(0.8, 0.05), Vec2.create(0.85, 0.05), Vec2.create(0.86, 0.0),
        ]),
        list: PD.ColorScale<ColorListName>('RedYellowBlue', ColorListOptions),
    }
    export type Params = typeof Params

    export const Utils: GeometryUtils<DirectVolume, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState,
        updateRenderableState
    }

    function createValues(directVolume: DirectVolume, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): DirectVolumeValues {
        const { gridTexture, gridTextureDim } = directVolume
        const { bboxSize, bboxMin, bboxMax, gridDimension, transform: gridTransform } = directVolume

        const { instanceCount, groupCount } = locationIt
        const color = createColors(locationIt, theme.color)
        const marker = createMarkers(instanceCount * groupCount)
        const overpaint = createEmptyOverpaint()
        const transparency = createEmptyTransparency()

        const counts = { drawCount: VolumeBox.indices.length, groupCount, instanceCount }

        const { boundingSphere, invariantBoundingSphere } = getBoundingSphere(gridDimension.ref.value, gridTransform.ref.value, transform.aTransform.ref.value, transform.instanceCount.ref.value)

        const controlPoints = getControlPointsFromVec2Array(props.controlPoints)
        const transferTex = createTransferFunctionTexture(controlPoints, props.list)

        const maxSteps = Math.ceil(Vec3.magnitude(gridDimension.ref.value)) * 2 * 5

        return {
            ...color,
            ...marker,
            ...overpaint,
            ...transparency,
            ...transform,
            ...BaseGeometry.createValues(props, counts),

            aPosition: ValueCell.create(VolumeBox.vertices as Float32Array),
            elements: ValueCell.create(VolumeBox.indices as Uint32Array),
            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),

            uIsoValue: ValueCell.create(props.isoValueNorm),
            uBboxMin: bboxMin,
            uBboxMax: bboxMax,
            uBboxSize: bboxSize,
            dMaxSteps: ValueCell.create(maxSteps),
            uTransform: gridTransform,
            uGridDim: gridDimension,
            dRenderMode: ValueCell.create(props.renderMode),
            tTransferTex: transferTex,

            dGridTexType: ValueCell.create(gridTexture.ref.value.depth > 0 ? '3d' : '2d'),
            uGridTexDim: gridTextureDim,
            tGridTex: gridTexture,
        }
    }

    function createValuesSimple(directVolume: DirectVolume, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform)
        const p = { ...PD.getDefaultValues(Params), ...props }
        return createValues(directVolume, s.transform, s.locationIterator, s.theme, p)
    }

    function updateValues(values: DirectVolumeValues, props: PD.Values<Params>) {
        ValueCell.updateIfChanged(values.uIsoValue, props.isoValueNorm)
        ValueCell.updateIfChanged(values.uAlpha, props.alpha)
        ValueCell.updateIfChanged(values.dUseFog, props.useFog)
        ValueCell.updateIfChanged(values.dRenderMode, props.renderMode)

        const controlPoints = getControlPointsFromVec2Array(props.controlPoints)
        createTransferFunctionTexture(controlPoints, props.list, values.tTransferTex)
    }

    function updateBoundingSphere(values: DirectVolumeValues, directVolume: DirectVolume) {
        const { boundingSphere, invariantBoundingSphere } = getBoundingSphere(values.uGridDim.ref.value, values.uTransform.ref.value, values.aTransform.ref.value, values.instanceCount.ref.value)
        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere)
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere)
        }
    }

    function createRenderableState(props: PD.Values<Params>): RenderableState {
        const state = BaseGeometry.createRenderableState(props)
        state.opaque = false
        return state
    }

    function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
        BaseGeometry.updateRenderableState(state, props)
        state.opaque = false
    }
}

//

const mTmp = Mat4.identity()
const mTmp2 = Mat4.identity()
const vHalfUnit = Vec3.create(0.5, 0.5, 0.5)
const tmpVertices = new Float32Array(VolumeBox.vertices.length)
function getBoundingSphere(gridDimension: Vec3, gridTransform: Mat4, transform: Float32Array, transformCount: number) {
    tmpVertices.set(VolumeBox.vertices)
    Mat4.fromTranslation(mTmp, vHalfUnit)
    Mat4.mul(mTmp, Mat4.fromScaling(mTmp2, gridDimension), mTmp)
    Mat4.mul(mTmp, gridTransform, mTmp)
    transformPositionArray(mTmp, tmpVertices, 0, tmpVertices.length / 3)
    return calculateBoundingSphere(tmpVertices, tmpVertices.length / 3, transform, transformCount)
}