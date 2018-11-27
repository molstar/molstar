/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, ElementIndex } from 'mol-model/structure';
import { GaussianDensity } from 'mol-math/geometry/gaussian-density';
import { Task, RuntimeContext } from 'mol-task';
import { DensityData } from 'mol-math/geometry';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { GaussianDensityTexture } from 'mol-math/geometry/gaussian-density/gpu';
import { Texture } from 'mol-gl/webgl/texture';
import { WebGLContext } from 'mol-gl/webgl/context';
import { PhysicalSizeTheme } from 'mol-theme/size/physical';

export const GaussianDensityParams = {
    resolution: PD.Numeric(1, { min: 0.1, max: 10, step: 0.1 }),
    radiusOffset: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }),
    smoothness: PD.Numeric(1.5, { min: 0.5, max: 2.5, step: 0.1 }),
    useGpu: PD.Boolean(false),
    ignoreCache: PD.Boolean(false),
}
export const DefaultGaussianDensityProps = PD.getDefaultValues(GaussianDensityParams)
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
    const sizeTheme = PhysicalSizeTheme({}, {})
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