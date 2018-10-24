/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData } from 'mol-model/volume'
import { RuntimeContext } from 'mol-task'
import { VolumeVisual, VolumeRepresentation } from '.';
import { DirectVolumeRenderObject, createDirectVolumeRenderObject } from 'mol-gl/render-object';
import { PickingId } from '../../geometry/picking';
import { MarkerAction } from '../../geometry/marker-data';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { createRenderableState, updateRenderableState, Geometry } from '../../geometry/geometry';
import { paramDefaultValues } from 'mol-view/parameter';
import { DirectVolume } from '../../geometry/direct-volume/direct-volume';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { Box3D } from 'mol-math/geometry';
import { Context } from 'mol-gl/webgl/context';
import { createTexture } from 'mol-gl/webgl/texture';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { NullLocation } from 'mol-model/location';
import { createIdentityTransform } from 'mol-geo/geometry/transform-data';

function getBoundingBox(gridDimension: Vec3, transform: Mat4) {
    const bbox = Box3D.empty()
    Box3D.add(bbox, gridDimension)
    Box3D.transform(bbox, bbox, transform)
    return bbox
}

// 2d volume texture

function getVolumeTexture2dLayout(dim: Vec3, maxTextureSize: number) {
    let width = 0
    let height = dim[1]
    let rows = 1
    let columns = dim[0]
    if (maxTextureSize < dim[0] * dim[2]) {
        columns =  Math.floor(maxTextureSize / dim[0])
        rows = Math.ceil(dim[2] / columns)
        width = columns * dim[0]
        height *= rows
    } else {
        width = dim[0] * dim[2]
    }
    width += columns // horizontal padding
    height += rows // vertical padding
    return { width, height, columns, rows }
}

function createVolumeTexture2d(volume: VolumeData, maxTextureSize: number) {
    const { data: tensor, dataStats: stats } = volume
    const { space, data } = tensor
    const dim = space.dimensions as Vec3
    const { get } = space
    const { width, height, columns, rows } = getVolumeTexture2dLayout(dim, maxTextureSize)

    const array = new Uint8Array(width * height * 4)
    const textureImage = { array, width, height }

    const [ xl, yl, zl ] = dim
    const xlp = xl + 1 // horizontal padding
    const ylp = yl + 1 // vertical padding

    function setTex(value: number, x: number, y: number, z: number) {
        const column = Math.floor(((z * xlp) % width) / xlp)
        const row = Math.floor((z * xlp) / width)
        const px = column * xlp + x
        const index = 4 * ((row * ylp * width) + (y * width) + px)
        array[index + 3] = ((value - stats.min) / (stats.max - stats.min)) * 255
    }

    console.log('dim', dim)
    console.log('layout', { width, height, columns, rows })

    for (let z = 0; z < zl; ++z) {
        for (let y = 0; y < yl; ++y) {
            for (let x = 0; x < xl; ++x) {
                setTex(get(data, x, y, z), x, y, z)
            }
        }
    }

    return textureImage
}

export function createDirectVolume2d(ctx: RuntimeContext, webgl: Context, volume: VolumeData, directVolume?: DirectVolume) {
    const gridDimension = volume.data.space.dimensions as Vec3
    const textureImage = createVolumeTexture2d(volume, webgl.maxTextureSize)
    // debugTexture(createImageData(textureImage.array, textureImage.width, textureImage.height), 1/3)
    const transform = VolumeData.getGridToCartesianTransform(volume)
    const bbox = getBoundingBox(gridDimension, transform)
    const dim = Vec3.create(gridDimension[0], gridDimension[1], gridDimension[2])
    dim[0] += 1 // horizontal padding
    dim[0] += 1 // vertical padding

    const texture = directVolume ? directVolume.gridTexture.ref.value : createTexture(webgl, 'image-uint8', 'rgba', 'ubyte', 'linear')
    texture.load(textureImage)

    return DirectVolume.create(bbox, dim, transform, texture, directVolume)
}

// 3d volume texture

function createVolumeTexture3d(volume: VolumeData) {
    const { data: tensor, dataStats: stats } = volume
    const { space, data } = tensor
    const [ width, height, depth ] = space.dimensions as Vec3
    const { get } = space

    const array = new Uint8Array(width * height * depth * 4)
    const textureVolume = { array, width, height, depth }

    let i = 0
    for (let z = 0; z < depth; ++z) {
        for (let y = 0; y < height; ++y) {
            for (let x = 0; x < width; ++x) {
                array[i + 3] = ((get(data, x, y, z) - stats.min) / (stats.max - stats.min)) * 255
                i += 4
            }
        }
    }

    return textureVolume
}

export function createDirectVolume3d(ctx: RuntimeContext, webgl: Context, volume: VolumeData, directVolume?: DirectVolume) {
    const gridDimension = volume.data.space.dimensions as Vec3
    const textureVolume = createVolumeTexture3d(volume)
    const transform = VolumeData.getGridToCartesianTransform(volume)
    const bbox = getBoundingBox(gridDimension, transform)

    const texture = directVolume ? directVolume.gridTexture.ref.value : createTexture(webgl, 'volume-uint8', 'rgba', 'ubyte', 'linear')
    texture.load(textureVolume)

    return DirectVolume.create(bbox, gridDimension, transform, texture, directVolume)
}

//

export const DirectVolumeParams = {
    ...Geometry.Params,
    ...DirectVolume.Params
}
export const DefaultDirectVolumeProps = paramDefaultValues(DirectVolumeParams)
export type DirectVolumeProps = typeof DefaultDirectVolumeProps

export function DirectVolumeVisual(): VolumeVisual<DirectVolumeProps> {
    let currentProps = DefaultDirectVolumeProps
    let renderObject: DirectVolumeRenderObject
    let currentVolume: VolumeData
    let directVolume: DirectVolume

    async function create(ctx: RuntimeContext, volume: VolumeData, props: Partial<DirectVolumeProps> = {}) {
        const { webgl } = props
        if (webgl === undefined) throw new Error('DirectVolumeVisual requires `webgl` in props')

        currentProps = { ...DefaultDirectVolumeProps, ...props }
        if (props.isoValueRelative) {
            // currentProps.isoValueAbsolute = VolumeIsoValue.calcAbsolute(currentVolume.dataStats, props.isoValueRelative)
        }

        const state = createRenderableState(currentProps)
        const locationIt = LocationIterator(1, 1, () => NullLocation)
        const transform = createIdentityTransform()

        if (webgl.isWebGL2) {
            directVolume = await createDirectVolume3d(ctx, webgl, volume, directVolume)
            const values = await DirectVolume.createValues(ctx, directVolume, transform, locationIt, currentProps)
            renderObject = createDirectVolumeRenderObject(values, state)
        } else {
            directVolume = await createDirectVolume2d(ctx, webgl, volume, directVolume)
            const values = await DirectVolume.createValues(ctx, directVolume, transform, locationIt, currentProps)
            renderObject = createDirectVolumeRenderObject(values, state)
        }
    }

    async function update(ctx: RuntimeContext, props: Partial<DirectVolumeProps> = {}) {
        const { webgl } = props
        if (webgl === undefined) throw new Error('DirectVolumeVisual requires `webgl` in props')

        const newProps = { ...currentProps, ...props }
        if (props.isoValueRelative) {
            // newProps.isoValueAbsolute = VolumeIsoValue.calcAbsolute(currentVolume.dataStats, props.isoValueRelative)
        }

        DirectVolume.updateValues(renderObject.values, newProps)
        updateRenderableState(renderObject.state, newProps)

        currentProps = newProps
    }

    return {
        get renderObject () { return renderObject },
        async createOrUpdate(ctx: RuntimeContext, props: Partial<DirectVolumeProps> = {}, volume?: VolumeData) {
            if (!volume && !currentVolume) {
                throw new Error('missing volume')
            } else if (volume && (!currentVolume || !renderObject)) {
                currentVolume = volume
                await create(ctx, volume, props)
            } else if (volume && volume !== currentVolume) {
                currentVolume = volume
                await create(ctx, volume, props)
            } else {
                await update(ctx, props)
            }

            currentProps = { ...DefaultDirectVolumeProps, ...props }
        },
        getLoci(pickingId: PickingId) {
            // TODO
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            // TODO
            return false
        },
        destroy() {
            // TODO
        }
    }
}

export function DirectVolumeRepresentation(): VolumeRepresentation<DirectVolumeProps> {
    let currentProps: DirectVolumeProps
    const volumeRepr = VolumeRepresentation(DirectVolumeVisual)
    return {
        label: 'Direct Volume',
        params: DirectVolumeParams,
        get renderObjects() {
            return [ ...volumeRepr.renderObjects ]
        },
        get props() {
            return { ...volumeRepr.props }
        },
        createOrUpdate: (props: Partial<DirectVolumeProps> = {}, volume?: VolumeData) => {
            currentProps = Object.assign({}, DefaultDirectVolumeProps, currentProps, props)
            return volumeRepr.createOrUpdate(currentProps, volume)
        },
        getLoci: (pickingId: PickingId) => {
            return volumeRepr.getLoci(pickingId)
        },
        mark: (loci: Loci, action: MarkerAction) => {
            return volumeRepr.mark(loci, action)
        },
        destroy() {
            volumeRepr.destroy()
        }
    }
}