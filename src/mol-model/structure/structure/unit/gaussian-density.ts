/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, ElementIndex } from 'mol-model/structure';
import { SizeTheme } from 'mol-view/theme/size';
import { GaussianDensity } from 'mol-math/geometry/gaussian-density';
import { Task, RuntimeContext } from 'mol-task';
import { DensityData } from 'mol-math/geometry';

export const DefaultGaussianDensityProps = {
    resolutionFactor: 6,
    radiusOffset: 0,
    smoothness: 1.5,
}
export type GaussianDensityProps = typeof DefaultGaussianDensityProps

function getConformation(unit: Unit) {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return unit.model.atomicConformation
        case Unit.Kind.Spheres: return unit.model.coarseConformation.spheres
        case Unit.Kind.Gaussians: return unit.model.coarseConformation.gaussians
    }
}

export function computeUnitGaussianDensity(unit: Unit, props: GaussianDensityProps) {
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

    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensity(ctx, position, unit.lookup3d.boundary.box, radius, { ...props, box: unit.lookup3d.boundary.box });
    });
}

export async function computeUnitGaussianDensityCached(unit: Unit, props: GaussianDensityProps, cache: Map<string, DensityData>, ctx?: RuntimeContext) {
    const key = `${props.radiusOffset}|${props.resolutionFactor}|${props.smoothness}`
    let density = cache.get(key)
    if (density) return density
    density = ctx ? await computeUnitGaussianDensity(unit, props).runInContext(ctx) : await computeUnitGaussianDensity(unit, props).run()
    cache.set(key, density)
    return density
}