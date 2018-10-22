/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext } from 'mol-task'
import { ValueCell } from 'mol-util'
import { Sphere3D, Box3D } from 'mol-math/geometry'
import { paramDefaultValues, RangeParam, SelectParam, TextParam } from 'mol-view/parameter';
import { DirectVolume2dValues, DirectVolumeBaseValues, DirectVolume3dValues } from 'mol-gl/renderable/direct-volume';
import { Vec3, Vec2, Mat4 } from 'mol-math/linear-algebra';
import { Box } from '../../primitive/box';
import { getControlPointsFromString, createTransferFunctionTexture } from './transfer-function';
import { Texture } from 'mol-gl/webgl/texture';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { TransformData } from '../transform-data';
import { createColors } from '../color-data';
import { createMarkers } from '../marker-data';
import { Geometry } from '../geometry';

const VolumeBox = Box()
const RenderModeOptions = [['isosurface', 'Isosurface'], ['volume', 'Volume']] as [string, string][]

interface DirectVolumeBase {
    readonly gridDimension: ValueCell<Vec3>,
    readonly bboxSize: ValueCell<Vec3>
    readonly bboxMin: ValueCell<Vec3>
    readonly bboxMax: ValueCell<Vec3>
    readonly transform: ValueCell<Mat4>

    /** Bounding sphere of the volume */
    boundingSphere?: Sphere3D
}

const BaseParams = {
    ...Geometry.Params,
    isoValueAbsolute: RangeParam('Iso Value Absolute', '', 0.22, -1, 1, 0.01),
    isoValueRelative: RangeParam('Iso Value Relative', '', 2, -10, 10, 0.1),
    renderMode: SelectParam('Render Mode', '', 'isosurface', RenderModeOptions),
    controlPoints: TextParam('Control Points', '', '0.19:0.1, 0.2:0.5, 0.21:0.1, 0.4:0.3'),
}
const DefaultBaseProps = paramDefaultValues(BaseParams)
type BaseProps = typeof DefaultBaseProps

async function createBaseValues(ctx: RuntimeContext, directVolume: DirectVolumeBase, transform: TransformData, locationIt: LocationIterator, props: BaseProps): Promise<DirectVolumeBaseValues> {
    const { instanceCount, groupCount } = locationIt
    const color = await createColors(ctx, locationIt, props)
    const marker = createMarkers(instanceCount * groupCount)

    const counts = { drawCount: VolumeBox.indices.length, groupCount, instanceCount }

    const { bboxSize, bboxMin, bboxMax, gridDimension, transform: gridTransform } = directVolume

    const controlPoints = getControlPointsFromString(props.controlPoints)
    const transferTex = createTransferFunctionTexture(controlPoints)

    const maxSteps = Math.ceil(Vec3.magnitude(gridDimension.ref.value)) * 2

    return {
        ...color,
        ...marker,
        ...transform,
        ...Geometry.createValues(props, counts),

        aPosition: ValueCell.create(VolumeBox.vertices as Float32Array),
        elements: ValueCell.create(VolumeBox.indices as Uint32Array),

        uIsoValue: ValueCell.create(props.isoValueAbsolute),
        uBboxMin: bboxMin,
        uBboxMax: bboxMax,
        uBboxSize: bboxSize,
        dMaxSteps: ValueCell.create(maxSteps),
        uTransform: gridTransform,
        uGridDim: gridDimension,
        dRenderMode: ValueCell.create(props.renderMode),
        tTransferTex: transferTex,
    }
}

function updateBaseValues(values: DirectVolumeBaseValues, props: BaseProps) {
    ValueCell.updateIfChanged(values.uIsoValue, props.isoValueAbsolute)
    ValueCell.updateIfChanged(values.uAlpha, props.alpha)
    ValueCell.updateIfChanged(values.dUseFog, props.useFog)
    ValueCell.updateIfChanged(values.dRenderMode, props.renderMode)

    const controlPoints = getControlPointsFromString(props.controlPoints)
    createTransferFunctionTexture(controlPoints, values.tTransferTex)
}

// 2d

export interface DirectVolume2d extends DirectVolumeBase {
    readonly kind: 'direct-volume-2d',
    readonly gridTexture: ValueCell<Texture>,
    readonly gridTextureDim: ValueCell<Vec2>,
}

export namespace DirectVolume2d {
    export function create(bbox: Box3D, gridDimension: Vec3, transform: Mat4, texture: Texture, directVolume?: DirectVolume2d): DirectVolume2d {
        if (directVolume) {
            ValueCell.update(directVolume.gridDimension, gridDimension)
            ValueCell.update(directVolume.gridTextureDim, Vec2.set(directVolume.gridTextureDim.ref.value, texture.width, texture.height))
            ValueCell.update(directVolume.bboxMin, bbox.min)
            ValueCell.update(directVolume.bboxMax, bbox.max)
            ValueCell.update(directVolume.bboxSize, Vec3.sub(directVolume.bboxSize.ref.value, bbox.max, bbox.min))
            ValueCell.update(directVolume.transform, transform)
            return directVolume
        } else {
            return {
                kind: 'direct-volume-2d' as 'direct-volume-2d',
                gridDimension: ValueCell.create(gridDimension),
                gridTexture: ValueCell.create(texture),
                gridTextureDim: ValueCell.create(Vec2.create(texture.width, texture.height)),
                bboxMin: ValueCell.create(bbox.min),
                bboxMax: ValueCell.create(bbox.max),
                bboxSize: ValueCell.create(Vec3.sub(Vec3.zero(), bbox.max, bbox.min)),
                transform: ValueCell.create(transform),
            }
        }
    }

    export function createEmpty(directVolume?: DirectVolume2d): DirectVolume2d {
        return {} as DirectVolume2d // TODO
    }

    export const Params = BaseParams
    export const DefaultProps = paramDefaultValues(Params)
    export type Props = typeof DefaultProps

    export async function createValues(ctx: RuntimeContext, directVolume: DirectVolume2d, transform: TransformData, locationIt: LocationIterator, props: Props): Promise<DirectVolume2dValues> {
        const { gridTexture, gridTextureDim } = directVolume

        return {
            ...await createBaseValues(ctx, directVolume, transform, locationIt, props),
            dGridTexType: ValueCell.create('2d'),
            uGridTexDim: gridTextureDim,
            tGridTex: gridTexture,
        }
    }

    export function updateValues(values: DirectVolume2dValues, props: Props) {
        updateBaseValues(values, props)
    }
}

// 3d

export interface DirectVolume3d extends DirectVolumeBase {
    readonly kind: 'direct-volume-3d',
    readonly gridTexture: ValueCell<Texture>,
}

export namespace DirectVolume3d {
    export function create(bbox: Box3D, gridDimension: Vec3, transform: Mat4, texture: Texture, directVolume?: DirectVolume3d): DirectVolume3d {
        if (directVolume) {
            ValueCell.update(directVolume.gridDimension, gridDimension)
            ValueCell.update(directVolume.bboxMin, bbox.min)
            ValueCell.update(directVolume.bboxMax, bbox.max)
            ValueCell.update(directVolume.bboxSize, Vec3.sub(directVolume.bboxSize.ref.value, bbox.max, bbox.min))
            ValueCell.update(directVolume.transform, transform)
            return directVolume
        } else {
            return {
                kind: 'direct-volume-3d' as 'direct-volume-3d',
                gridDimension: ValueCell.create(gridDimension),
                gridTexture: ValueCell.create(texture),
                bboxMin: ValueCell.create(bbox.min),
                bboxMax: ValueCell.create(bbox.max),
                bboxSize: ValueCell.create(Vec3.sub(Vec3.zero(), bbox.max, bbox.min)),
                transform: ValueCell.create(transform),
            }
        }
    }

    export function createEmpty(directVolume?: DirectVolume3d): DirectVolume3d {
        return {} as DirectVolume3d // TODO
    }

    export const Params = BaseParams
    export const DefaultProps = paramDefaultValues(Params)
    export type Props = typeof DefaultProps

    export async function createValues(ctx: RuntimeContext, directVolume: DirectVolume3d, transform: TransformData, locationIt: LocationIterator, props: Props): Promise<DirectVolume3dValues> {
        const { gridTexture } = directVolume

        return {
            ...await createBaseValues(ctx, directVolume, transform, locationIt, props),
            dGridTexType: ValueCell.create('3d'),
            tGridTex: gridTexture,
        }
    }

    export function updateValues(values: DirectVolume3dValues, props: Props) {
        updateBaseValues(values, props)
    }
}