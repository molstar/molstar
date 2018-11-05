/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, ElementIndex } from 'mol-model/structure';
import { SizeTheme } from 'mol-theme/size';
import { GaussianDensity } from 'mol-math/geometry/gaussian-density';
import { Task, RuntimeContext } from 'mol-task';
import { DensityData } from 'mol-math/geometry';
import { NumberParam, paramDefaultValues, BooleanParam } from 'mol-util/parameter';
import { GaussianDensityTexture } from 'mol-math/geometry/gaussian-density/gpu';
import { Texture } from 'mol-gl/webgl/texture';
import { WebGLContext } from 'mol-gl/webgl/context';

export const GaussianDensityParams = {
    resolution: NumberParam('Resolution', '', 1, 0.1, 10, 0.1),
    radiusOffset: NumberParam('Radius Offset', '', 0, 0, 10, 0.1),
    smoothness: NumberParam('Smoothness', '', 1.5, 0.5, 2.5, 0.1),
    useGpu: BooleanParam('Use GPU', '', true),
    ignoreCache: BooleanParam('Ignore Cache', '', false),
}
export const DefaultGaussianDensityProps = paramDefaultValues(GaussianDensityParams)
export type GaussianDensityProps = typeof DefaultGaussianDensityProps

function getConformation(unit: Unit) {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return unit.model.atomicConformation
        case Unit.Kind.Spheres: return unit.model.coarseConformation.spheres
        case Unit.Kind.Gaussians: return unit.model.coarseConformation.gaussians
    }
}

function getConformationAndRadius(unit: Unit) {
    const conformation = getConformation(unit)
    const { elements } = unit
    const position = {
        indices: elements,
        x: conformation.x,
        y: conformation.y,
        z: conformation.z
    }

    const l = StructureElement.create(unit)
    const sizeTheme = SizeTheme({ name: 'physical' })
    const radius = (index: number) => {
        l.element = index as ElementIndex
        return sizeTheme.size(l)
    }

    return { position, radius }
}

export function computeUnitGaussianDensity(unit: Unit, props: GaussianDensityProps, webgl?: WebGLContext) {
    const { position, radius } = getConformationAndRadius(unit)
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensity(ctx, position, unit.lookup3d.boundary.box, radius, props, webgl);
    });
}

export function computeUnitGaussianDensityTexture(unit: Unit, props: GaussianDensityProps, webgl: WebGLContext, texture?: Texture) {
    const { position, radius } = getConformationAndRadius(unit)
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensityTexture(ctx, webgl, position, unit.lookup3d.boundary.box, radius, props, texture);
    });
}

export async function computeUnitGaussianDensityCached(unit: Unit, props: GaussianDensityProps, cache: Map<string, DensityData>, ctx: RuntimeContext, webgl?: WebGLContext) {
    const key = `${props.radiusOffset}|${props.resolution}|${props.smoothness}`
    let density = cache.get(key)
    if (density && !props.ignoreCache) return density
    density = await computeUnitGaussianDensity(unit, props, webgl).runInContext(ctx)
    if (!props.ignoreCache) cache.set(key, density)
    return density
}