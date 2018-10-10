/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext } from 'mol-task'
import { ValueCell } from 'mol-util'
import { Sphere3D } from 'mol-math/geometry'
import { paramDefaultValues, RangeParam, BooleanParam, SelectParam, TextParam } from 'mol-view/parameter';
import { DirectVolume2dValues, DirectVolumeBaseValues, DirectVolume3dValues } from 'mol-gl/renderable/direct-volume';
import { TextureImage, TextureVolume } from 'mol-gl/renderable/util';
import { Vec3, Vec2, Mat4 } from 'mol-math/linear-algebra';
import { Box } from '../../primitive/box';
import { getControlPointsFromString, createTransferFunctionTexture } from './transfer-function';

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
    alpha: RangeParam('Opacity', '', 1, 0, 1, 0.01),
    visible: BooleanParam('Visible', '', true),
    depthMask: BooleanParam('Depth Mask', '', true),
    useFog: BooleanParam('Use Fog', '', false),
    isoValueAbsolute: RangeParam('Iso Value Absolute', '', 0.22, -1, 1, 0.01),
    isoValueRelative: RangeParam('Iso Value Relative', '', 2, -10, 10, 0.1),
    renderMode: SelectParam('Render Mode', '', 'volume', RenderModeOptions),
    controlPoints: TextParam('Control Points', '', '0.19:0.1, 0.2:0.5, 0.21:0.1, 0.4:0.3'),
}
const DefaultBaseProps = paramDefaultValues(BaseParams)
type BaseProps = typeof DefaultBaseProps

async function createBaseValues(ctx: RuntimeContext, directVolume: DirectVolumeBase, props: BaseProps): Promise<DirectVolumeBaseValues> {
    const { bboxSize, bboxMin, bboxMax, gridDimension, transform } = directVolume

    const controlPoints = getControlPointsFromString(props.controlPoints)
    const transferTex = createTransferFunctionTexture(controlPoints)

    const maxSteps = Math.round(Vec3.magnitude(bboxSize.ref.value))
    console.log('maxSteps', maxSteps)

    return {
        drawCount: ValueCell.create(VolumeBox.indices.length),
        instanceCount: ValueCell.create(1),

        aPosition: ValueCell.create(VolumeBox.vertices as Float32Array),
        elements: ValueCell.create(VolumeBox.indices as Uint32Array),

        uAlpha: ValueCell.create(props.alpha),
        dUseFog: ValueCell.create(props.useFog),

        uIsoValue: ValueCell.create(props.isoValueAbsolute),
        uBboxMin: bboxMin,
        uBboxMax: bboxMax,
        uBboxSize: bboxSize,
        dMaxSteps: ValueCell.create(maxSteps),
        uTransform: transform,
        uGridDim: gridDimension,
        dRenderMode: ValueCell.create(props.renderMode),
        tTransferTex: transferTex,
    }
}

function updateBaseValues(values: DirectVolumeBaseValues, props: BaseProps) {
    console.log('DirectVolumeBaseValues', props, values)
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
    readonly gridTexture: ValueCell<TextureImage<any>>,
    readonly gridTextureDim: ValueCell<Vec2>,
}

export namespace DirectVolume2d {
    export function createEmpty(directVolume?: DirectVolume2d): DirectVolume2d {
        return {} as DirectVolume2d // TODO
    }

    export const Params = BaseParams
    export const DefaultProps = paramDefaultValues(Params)
    export type Props = typeof DefaultProps

    export async function createValues(ctx: RuntimeContext, directVolume: DirectVolume2d, props: Props): Promise<DirectVolume2dValues> {
        const { gridTexture, gridTextureDim } = directVolume

        return {
            ...await createBaseValues(ctx, directVolume, props),
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
    readonly gridTexture: ValueCell<TextureVolume<any>>,
}

export namespace DirectVolume3d {
    export function createEmpty(directVolume?: DirectVolume3d): DirectVolume3d {
        return {} as DirectVolume3d // TODO
    }

    export const Params = BaseParams
    export const DefaultProps = paramDefaultValues(Params)
    export type Props = typeof DefaultProps

    export async function createValues(ctx: RuntimeContext, directVolume: DirectVolume3d, props: Props): Promise<DirectVolume3dValues> {
        const { gridTexture } = directVolume

        return {
            ...await createBaseValues(ctx, directVolume, props),
            dGridTexType: ValueCell.create('3d'),
            tGridTex: gridTexture,
        }
    }

    export function updateValues(values: DirectVolume3dValues, props: Props) {
        updateBaseValues(values, props)
    }
}