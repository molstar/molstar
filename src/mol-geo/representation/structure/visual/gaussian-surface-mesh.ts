/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure, StructureElement, ElementIndex } from 'mol-model/structure';
import { UnitsVisual, VisualUpdateState } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../geometry/mesh/mesh';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { computeMarchingCubesMesh } from '../../../util/marching-cubes/algorithm';
import { GaussianDensityProps, GaussianDensityParams } from 'mol-model/structure/structure/unit/gaussian-density';
import { paramDefaultValues } from 'mol-view/parameter';
import { SizeTheme } from 'mol-view/theme/size';
import { OrderedSet } from 'mol-data/int';

async function createGaussianSurfaceMesh(ctx: RuntimeContext, unit: Unit, structure: Structure, props: GaussianDensityProps, mesh?: Mesh): Promise<Mesh> {
    const { smoothness, radiusOffset } = props
    const { transform, field, idField } = await unit.computeGaussianDensity(props, ctx)

    const params = {
        isoLevel: Math.exp(-smoothness),
        scalarField: field,
        idField
    }
    const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx)

    Mesh.transformImmediate(surface, transform)

    if (props.useGpu) {
        console.time('find max element radius')
        const { elements } = unit
        const n = OrderedSet.size(elements)
        const l = StructureElement.create(unit)
        const sizeTheme = SizeTheme({ name: 'physical' })
        const radius = (index: number) => {
            l.element = index as ElementIndex
            return sizeTheme.size(l)
        }
        let maxRadius = 0
        for (let i = 0; i < n; ++i) {
            const r = radius(OrderedSet.getAt(elements, i)) + radiusOffset
            if (maxRadius < r) maxRadius = r
        }
        console.timeEnd('find max element radius')

        console.time('find closest element for vertices')
        const { lookup3d } = unit

        const { vertexCount, vertexBuffer, groupBuffer } = surface
        const vertices = vertexBuffer.ref.value
        const groups = groupBuffer.ref.value
        for (let i = 0; i < vertexCount; ++i) {
            const r = lookup3d.find(vertices[i * 3], vertices[i * 3 + 1], vertices[i * 3 + 2], maxRadius * 2)
            let minDsq = Infinity
            let group = 0
            for (let j = 0, jl = r.count; j < jl; ++j) {
                const dSq = r.squaredDistances[j]
                if (dSq < minDsq) {
                    minDsq = dSq
                    group = r.indices[j]
                }
            }
            groups[i] = group
        }
        console.timeEnd('find closest element for vertices')
    }

    Mesh.computeNormalsImmediate(surface)
    Mesh.uniformTriangleGroup(surface)

    return surface;
}

export const GaussianSurfaceParams = {
    ...UnitsMeshParams,
    ...GaussianDensityParams,
}
export const DefaultGaussianSurfaceProps = paramDefaultValues(GaussianSurfaceParams)
export type GaussianSurfaceProps = typeof DefaultGaussianSurfaceProps

export function GaussianSurfaceVisual(): UnitsVisual<GaussianSurfaceProps> {
    return UnitsMeshVisual<GaussianSurfaceProps>({
        defaultProps: DefaultGaussianSurfaceProps,
        createGeometry: createGaussianSurfaceMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: GaussianSurfaceProps, currentProps: GaussianSurfaceProps) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true
            if (newProps.useGpu !== currentProps.useGpu) state.createGeometry = true
            if (newProps.ignoreCache !== currentProps.ignoreCache) state.createGeometry = true
        }
    })
}