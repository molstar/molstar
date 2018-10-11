/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData } from 'mol-model/volume'
import { RuntimeContext } from 'mol-task'
import { VolumeVisual, VolumeRepresentation } from '.';
import { DirectVolume2dRenderObject, createDirectVolume2dRenderObject, DirectVolume3dRenderObject, createDirectVolume3dRenderObject } from 'mol-gl/render-object';
import { PickingId } from '../../geometry/picking';
import { MarkerAction } from '../../geometry/marker-data';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { createRenderableState, updateRenderableState, Geometry } from '../../geometry/geometry';
import { paramDefaultValues } from 'mol-view/parameter';
import { ValueCell } from 'mol-util';
import { DirectVolume2d, DirectVolume3d } from '../../geometry/direct-volume/direct-volume';
import { Vec2, Vec3 } from 'mol-math/linear-algebra';
import { Box3D } from 'mol-math/geometry';
import { createImageData } from 'mol-gl/webgl/context';
import { debugTexture } from 'mol-gl/util';
import { DirectVolume3dValues, DirectVolume2dValues } from 'mol-gl/renderable/direct-volume';

// 2d volume texture

function getVolumeTexture2dLayout(dim: Vec3, maxTextureSize = 4096) {
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

function createVolumeTexture2d(volume: VolumeData) {
    const { data: tensor, dataStats: stats } = volume
    const { space, data } = tensor
    const dim = space.dimensions as Vec3
    const { get } = space
    const { width, height, columns, rows } = getVolumeTexture2dLayout(dim)

    const array = new Uint8Array(width * height * 4)
    const textureImage = { array, width, height }

    const [ xl, yl, zl ] = dim
    const xlp = xl + 1
    const ylp = yl + 1

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

export function createDirectVolume2d(ctx: RuntimeContext, volume: VolumeData, directVolume?: DirectVolume2d) {
    const gridDimension = volume.data.space.dimensions as Vec3
    // const textureImage = createTextureImage(1, 4)
    const textureImage = createVolumeTexture2d(volume)
    const transform = VolumeData.getGridToCartesianTransform(volume)

    console.log('textureImage', textureImage)
    debugTexture(createImageData(textureImage.array, textureImage.width, textureImage.height), 1/3)

    const bbox = Box3D.empty()
    Box3D.add(bbox, gridDimension)
    Box3D.transform(bbox, bbox, transform)

    const dim = Vec3.create(gridDimension[0] + 1, gridDimension[1] + 1, gridDimension[2])

    if (directVolume) {
        ValueCell.update(directVolume.gridDimension, dim)
        ValueCell.update(directVolume.gridTexture, textureImage)
        ValueCell.update(directVolume.gridTextureDim, Vec2.set(directVolume.gridTextureDim.ref.value, textureImage.width, textureImage.height))
        ValueCell.update(directVolume.bboxMin, bbox.min)
        ValueCell.update(directVolume.bboxMax, bbox.max)
        ValueCell.update(directVolume.bboxSize, Vec3.sub(directVolume.bboxSize.ref.value, bbox.max, bbox.min))
        ValueCell.update(directVolume.transform, transform)
    } else {
        directVolume = {
            kind: 'direct-volume-2d' as 'direct-volume-2d',
            gridDimension: ValueCell.create(dim),
            gridTexture: ValueCell.create(textureImage),
            gridTextureDim: ValueCell.create(Vec2.create(textureImage.width, textureImage.height)),
            bboxMin: ValueCell.create(bbox.min),
            bboxMax: ValueCell.create(bbox.max),
            bboxSize: ValueCell.create(Vec3.sub(Vec3.zero(), bbox.max, bbox.min)),
            transform: ValueCell.create(transform),
        }
    }

    console.log('gridDimension', dim)
    console.log('gridTextureDim', textureImage.width, textureImage.height)
    console.log('boundingBox', bbox)
    console.log('transform', transform)

    return directVolume;
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

export function createDirectVolume3d(ctx: RuntimeContext, volume: VolumeData, directVolume?: DirectVolume3d) {
    const gridDimension = volume.data.space.dimensions as Vec3
    const textureVolume = createVolumeTexture3d(volume)
    const transform = VolumeData.getGridToCartesianTransform(volume)

    console.log('textureVolume', textureVolume)

    const bbox = Box3D.empty()
    Box3D.add(bbox, gridDimension)
    Box3D.transform(bbox, bbox, transform)

    if (directVolume) {
        ValueCell.update(directVolume.gridDimension, gridDimension)
        ValueCell.update(directVolume.gridTexture, textureVolume)
        ValueCell.update(directVolume.bboxMin, bbox.min)
        ValueCell.update(directVolume.bboxMax, bbox.max)
        ValueCell.update(directVolume.bboxSize, Vec3.sub(directVolume.bboxSize.ref.value, bbox.max, bbox.min))
        ValueCell.update(directVolume.transform, transform)
    } else {
        directVolume = {
            kind: 'direct-volume-3d' as 'direct-volume-3d',
            gridDimension: ValueCell.create(gridDimension),
            gridTexture: ValueCell.create(textureVolume),
            bboxMin: ValueCell.create(bbox.min),
            bboxMax: ValueCell.create(bbox.max),
            bboxSize: ValueCell.create(Vec3.sub(Vec3.zero(), bbox.max, bbox.min)),
            transform: ValueCell.create(transform),
        }
    }

    return directVolume;
}

//

function hasWebGL2() {
    return true
}

export const DirectVolumeParams = {
    ...Geometry.Params,
    ...DirectVolume2d.Params
}
export const DefaultDirectVolumeProps = paramDefaultValues(DirectVolumeParams)
export type DirectVolumeProps = typeof DefaultDirectVolumeProps

export function DirectVolumeVisual(): VolumeVisual<DirectVolumeProps> {
    let currentProps = DefaultDirectVolumeProps
    let renderObject: DirectVolume2dRenderObject | DirectVolume3dRenderObject
    let currentVolume: VolumeData
    let directVolume: DirectVolume2d | DirectVolume3d

    async function create(ctx: RuntimeContext, volume: VolumeData, props: Partial<DirectVolumeProps> = {}) {
        currentProps = { ...DefaultDirectVolumeProps, ...props }
        if (props.isoValueRelative) {
            // currentProps.isoValueAbsolute = VolumeIsoValue.calcAbsolute(currentVolume.dataStats, props.isoValueRelative)
        }

        const state = createRenderableState(currentProps)

        if (hasWebGL2()) {
            console.log('createing 3d volume')
            directVolume = await createDirectVolume3d(ctx, volume, directVolume as DirectVolume3d)
            const values = await DirectVolume3d.createValues(ctx, directVolume as DirectVolume3d, currentProps)
            renderObject = createDirectVolume3dRenderObject(values, state)
        } else {
            directVolume = await createDirectVolume2d(ctx, volume, directVolume as DirectVolume2d)
            const values = await DirectVolume2d.createValues(ctx, directVolume as DirectVolume2d, currentProps)
            renderObject = createDirectVolume2dRenderObject(values, state)
        }
    }

    async function update(ctx: RuntimeContext, props: Partial<DirectVolumeProps> = {}) {
        console.log('props', props)
        const newProps = { ...currentProps, ...props }
        if (props.isoValueRelative) {
            // newProps.isoValueAbsolute = VolumeIsoValue.calcAbsolute(currentVolume.dataStats, props.isoValueRelative)
        }

        if (hasWebGL2()) {
            DirectVolume3d.updateValues(renderObject.values as DirectVolume3dValues, newProps)
        } else {
            DirectVolume2d.updateValues(renderObject.values as DirectVolume2dValues, newProps)
        }
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