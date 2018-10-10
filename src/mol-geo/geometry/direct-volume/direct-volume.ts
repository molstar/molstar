/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext } from 'mol-task'
import { ValueCell } from 'mol-util'
import { Sphere3D } from 'mol-math/geometry'
import { paramDefaultValues, RangeParam, BooleanParam, SelectParam, TextParam } from 'mol-view/parameter';
import { DirectVolume2dValues } from 'mol-gl/renderable/direct-volume';
import { TextureImage } from 'mol-gl/renderable/util';
import { Vec3, Vec2, Mat4 } from 'mol-math/linear-algebra';
import { Box } from '../../primitive/box';
import { getControlPointsFromString, createTransferFunctionTexture } from './transfer-function';

export interface DirectVolume {
    readonly kind: 'direct-volume',

    readonly gridDimension: ValueCell<Vec3>,
    readonly gridTexture: ValueCell<TextureImage<any>>,
    readonly gridTextureDim: ValueCell<Vec2>,
    readonly bboxSize: ValueCell<Vec3>
    readonly bboxMin: ValueCell<Vec3>
    readonly bboxMax: ValueCell<Vec3>
    readonly transform: ValueCell<Mat4>

    /** Bounding sphere of the volume */
    boundingSphere?: Sphere3D
}

const VolumeBox = Box()

const RenderModeOptions = [['isosurface', 'Isosurface'], ['volume', 'Volume']] as [string, string][]

export namespace DirectVolume {
    export function createEmpty(directVolume?: DirectVolume): DirectVolume {
        // TODO
        return {

        } as DirectVolume
    }

    //

    export const Params = {
        alpha: RangeParam('Opacity', '', 1, 0, 1, 0.01),
        visible: BooleanParam('Visible', '', true),
        depthMask: BooleanParam('Depth Mask', '', true),
        useFog: BooleanParam('Use Fog', '', false),
        isoValueAbsolute: RangeParam('Iso Value Absolute', '', 0.22, -1, 1, 0.01),
        isoValueRelative: RangeParam('Iso Value Relative', '', 2, -10, 10, 0.1),
        renderMode: SelectParam('Render Mode', '', 'volume', RenderModeOptions),
        controlPoints: TextParam('Control Points', '', '0.19:0.1, 0.2:0.5, 0.21:0.1, 0.4:0.3'),
    }
    export const DefaultProps = paramDefaultValues(Params)
    export type Props = typeof DefaultProps

    export async function createValues(ctx: RuntimeContext, directVolume: DirectVolume, props: Props): Promise<DirectVolume2dValues> {
        const { bboxSize, bboxMin, bboxMax, gridDimension, gridTexture, gridTextureDim, transform } = directVolume

        const controlPoints = getControlPointsFromString(props.controlPoints)
        const transferTex = createTransferFunctionTexture(controlPoints)

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
            uTransform: transform,
            uGridDim: gridDimension,
            uGridTexDim: gridTextureDim,
            tGridTex: gridTexture,
            dRenderMode: ValueCell.create(props.renderMode),
            tTransferTex: transferTex,
        }
    }

    export function updateValues(values: DirectVolume2dValues, props: Props) {
        console.log('DirectVolumeValues', props, values)
        ValueCell.updateIfChanged(values.uIsoValue, props.isoValueAbsolute)
        ValueCell.updateIfChanged(values.uAlpha, props.alpha)
        ValueCell.updateIfChanged(values.dUseFog, props.useFog)
        ValueCell.updateIfChanged(values.dRenderMode, props.renderMode)

        const controlPoints = getControlPointsFromString(props.controlPoints)
        createTransferFunctionTexture(controlPoints, values.tTransferTex)
    }
}