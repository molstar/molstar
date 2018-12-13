/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, StructureElement } from 'mol-model/structure';
import { VisualUpdateState } from '../../util';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { DirectVolume } from 'mol-geo/geometry/direct-volume/direct-volume';
import { VisualContext } from 'mol-repr/visual';
import { Theme } from 'mol-theme/theme';
import { ComplexVisual, ComplexDirectVolumeVisual, ComplexDirectVolumeParams } from '../complex-visual';
import { EmptyLoci } from 'mol-model/loci';
import { NullLocation } from 'mol-model/location';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { WebGLContext } from 'mol-gl/webgl/context';
import { Texture } from 'mol-gl/webgl/texture';
import { GaussianDensityTexture } from 'mol-math/geometry/gaussian-density/gpu';
import { Task } from 'mol-task';
import { OrderedSet } from 'mol-data/int';
import { PhysicalSizeTheme } from 'mol-theme/size/physical';

function getConformationAndRadius(structure: Structure) {
    const n = structure.elementCount

    const xs = new Float32Array(n)
    const ys = new Float32Array(n)
    const zs = new Float32Array(n)
    const rs = new Float32Array(n)

    const l = StructureElement.create()
    const sizeTheme = PhysicalSizeTheme({}, {})

    let m = 0
    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const unit = structure.units[i]
        const { elements } = unit
        const { x, y, z } = unit.conformation
        l.unit = unit
        for (let j = 0, jl = elements.length; j < jl; ++j) {
            const eI = elements[j]
            xs[m + j] = x(eI)
            ys[m + j] = y(eI)
            zs[m + j] = z(eI)
            l.element = eI
            rs[m + j] = sizeTheme.size(l)
        }
        m += elements.length
    }

    const position = { indices: OrderedSet.ofRange(0, n), x: xs, y: ys, z: zs }
    const radius = (index: number) => rs[index]

    return { position, radius }
}

// TODO calculate by combining unit volumes
function computeStructureGaussianDensityTexture(structure: Structure, props: PD.Values<GaussianDensityVolumeParams>, webgl: WebGLContext, texture?: Texture) {
    const { position, radius } = getConformationAndRadius(structure)
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensityTexture(ctx, webgl, position, structure.lookup3d.boundary.box, radius, props, texture);
    });
}

async function createGaussianDensityVolume(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<GaussianDensityVolumeParams>, directVolume?: DirectVolume): Promise<DirectVolume> {
    const { runtime, webgl } = ctx
    if (webgl === undefined) throw new Error('createGaussianDensityVolume requires `webgl` object in VisualContext')

    const p = { ...props, useGpu: true }
    const oldTexture = directVolume ? directVolume.gridTexture.ref.value : undefined
    const densityTextureData = await computeStructureGaussianDensityTexture(structure, p, webgl, oldTexture).runInContext(runtime)
    const { transform, texture, bbox, gridDimension } = densityTextureData

    return DirectVolume.create(bbox, gridDimension, transform, texture, directVolume)
}

export const GaussianDensityVolumeParams = {
    ...ComplexDirectVolumeParams,
    resolution: PD.Numeric(1, { min: 0.1, max: 10, step: 0.1 }),
    radiusOffset: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }),
    smoothness: PD.Numeric(1.5, { min: 0.5, max: 2.5, step: 0.1 }),
}
export type GaussianDensityVolumeParams = typeof GaussianDensityVolumeParams

export function GaussianDensityVolumeVisual(): ComplexVisual<GaussianDensityVolumeParams> {
    return ComplexDirectVolumeVisual<GaussianDensityVolumeParams>({
        defaultProps: PD.getDefaultValues(GaussianDensityVolumeParams),
        createGeometry: createGaussianDensityVolume,
        createLocationIterator: (structure: Structure) => LocationIterator(structure.elementCount, 1, () => NullLocation),
        getLoci: () => EmptyLoci, // TODO
        mark: () => false, // TODO
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<GaussianDensityVolumeParams>, currentProps: PD.Values<GaussianDensityVolumeParams>) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true
            if (newProps.smoothness !== currentProps.smoothness) {
                state.createGeometry = true
                newProps.isoValue = Math.exp(-newProps.smoothness)
            }
        }
    })
}