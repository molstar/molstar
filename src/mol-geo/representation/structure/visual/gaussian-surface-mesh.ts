/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure, StructureElement } from 'mol-model/structure';
import { UnitsVisual, MeshUpdateState } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../mesh/mesh';
import { UnitsMeshVisual, DefaultUnitsMeshProps } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { computeMarchingCubes } from '../../../util/marching-cubes/algorithm';
import { Tensor, Vec3, Mat4 } from 'mol-math/linear-algebra';
import { Box3D } from 'mol-math/geometry';
import { SizeThemeProps, SizeTheme } from 'mol-view/theme/size';
import { Color } from 'mol-util/color';

export interface GaussianSurfaceMeshProps {
    sizeTheme: SizeThemeProps

    resolutionFactor: number
    probeRadius: number
    isoValue: number
}

function getDelta(box: Box3D, resolutionFactor: number) {
    const extent = Vec3.sub(Vec3.zero(), box.max, box.min)

    const n = Math.pow(Math.pow(2, resolutionFactor), 3)
    const f = (extent[0] * extent[1] * extent[2]) / n
    const s = Math.pow(f, 1 / 3)
    const size = Vec3.zero()
    // Vec3.scale(size, extent, s)
    Vec3.ceil(size, Vec3.scale(size, extent, s))
    const delta = Vec3.div(Vec3.zero(), extent, size)
    return delta
}

async function createGaussianSurfaceMesh(ctx: RuntimeContext, unit: Unit, structure: Structure, props: GaussianSurfaceMeshProps, mesh?: Mesh): Promise<Mesh> {
    const { resolutionFactor, probeRadius, isoValue } = props

    const { elements } = unit;
    const elementCount = elements.length;
    const sizeTheme = SizeTheme(props.sizeTheme)

    const v = Vec3.zero()
    const p = Vec3.zero()
    const pos = unit.conformation.invariantPosition
    const l = StructureElement.create()
    l.unit = unit

    const pad = (probeRadius + 2) * 2 // TODO calculate max radius
    const box = unit.lookup3d.boundary.box
    const expandedBox = Box3D.expand(Box3D.empty(), box, Vec3.create(pad, pad, pad));
    const extent = Vec3.sub(Vec3.zero(), expandedBox.max, expandedBox.min)
    const min = expandedBox.min

    // const n = Math.pow(128, 3)
    // const f = (extent[0] * extent[1] * extent[2]) / n
    // const s = Math.pow(f, 1 / 3)
    // const size = Vec3.zero()
    // // Vec3.scale(size, extent, s)
    // Vec3.ceil(size, Vec3.scale(size, extent, s))
    // const delta = Vec3.div(Vec3.zero(), extent, size)

    // console.log('extent', extent)
    // console.log('n', n)
    // console.log('f', f)
    // console.log('s', s)
    // console.log('size', size)
    // console.log('delta', delta)
    const delta = getDelta(Box3D.expand(Box3D.empty(), structure.boundary.box, Vec3.create(pad, pad, pad)), resolutionFactor)
    const dim = Vec3.zero()
    Vec3.ceil(dim, Vec3.mul(dim, extent, delta))
    // console.log('dim', dim, dim[0] * dim[1] * dim[2])

    const space = Tensor.Space(dim, [0, 1, 2], Float32Array)
    const data = space.create()
    const field = Tensor.create(space, data)

    const c = Vec3.zero()

    const radius2 = Vec3.scale(Vec3.zero(), delta, 2)
    const alpha = 1.5 // smoothness

    console.time('surface density')
    for (let i = 0; i < elementCount; i++) {
        l.element = elements[i]
        pos(elements[i], v)

        Vec3.sub(v, v, min)
        Vec3.mul(c, v, delta)

        const radius = sizeTheme.size(l) // + probeRadius

        const minX = Math.floor(c[0] - radius2[0]) - 1
        const minY = Math.floor(c[1] - radius2[1]) - 1
        const minZ = Math.floor(c[2] - radius2[2]) - 1
        const maxX = Math.ceil(c[0] + radius2[0]) + 1
        const maxY = Math.ceil(c[1] + radius2[1]) + 1
        const maxZ = Math.ceil(c[2] + radius2[2]) + 1

        for (let x = minX; x < maxX; ++x) {
            for (let y = minY; y < maxY; ++y) {
                for (let z = minZ; z < maxZ; ++z) {
                    Vec3.set(p, x, y, z)
                    Vec3.div(p, p, delta)
                    const dist = Vec3.distance(p, v)
                    const density = Math.exp(-alpha * Math.pow((dist) / radius, 2))
                    space.set(data, x, y, z, space.get(data, x, y, z) + density)
                }
            }
        }

        if (i % 100 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Gaussian surface', current: i, max: elementCount });
        }
    }
    console.timeEnd('surface density')

    // console.log('data', data)

    console.time('surface mc')
    const surface = await computeMarchingCubes({
        isoLevel: isoValue * Math.exp(-alpha),
        scalarField: field,
        oldSurface: mesh
    }).runAsChild(ctx)
    console.timeEnd('surface mc')

    const t = Mat4.identity()
    Mat4.fromScaling(t, Vec3.inverse(Vec3.zero(), delta))
    Mat4.setTranslation(t, expandedBox.min)

    Mesh.transformImmediate(surface, t)
    Mesh.computeNormalsImmediate(surface)

    return surface;
}

export const DefaultGaussianSurfaceProps = {
    ...DefaultUnitsMeshProps,
    linearSegments: 8,
    radialSegments: 12,
    aspectRatio: 5,
    arrowFactor: 1.5,

    flipSided: true,
    // flatShaded: true,
    alpha: 0.7,
    colorTheme: { name: 'uniform' as 'uniform', value: Color(0xDDDDDD) },

    resolutionFactor: 7,
    probeRadius: 1.4,
    isoValue: 1,
}
export type GaussianSurfaceProps = typeof DefaultGaussianSurfaceProps

export function GaussianSurfaceVisual(): UnitsVisual<GaussianSurfaceProps> {
    return UnitsMeshVisual<GaussianSurfaceProps>({
        defaultProps: DefaultGaussianSurfaceProps,
        createMesh: createGaussianSurfaceMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: MeshUpdateState, newProps: GaussianSurfaceProps, currentProps: GaussianSurfaceProps) => {
            if (newProps.resolutionFactor !== currentProps.resolutionFactor) state.createMesh = true
            if (newProps.probeRadius !== currentProps.probeRadius) state.createMesh = true
            if (newProps.isoValue !== currentProps.isoValue) state.createMesh = true
        }
    })
}