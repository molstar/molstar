/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual, MeshUpdateState } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../mesh/mesh';
import { UnitsMeshVisual, DefaultUnitsMeshProps } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { computeMarchingCubes } from '../../../util/marching-cubes/algorithm';
import { Tensor, Vec3, Mat4 } from 'mol-math/linear-algebra';
import { Box3D } from 'mol-math/geometry';
import { ValueCell } from 'mol-util';
import { smoothstep } from 'mol-math/interpolate';

export interface GaussianSurfaceMeshProps {
    
}

function getDelta(box: Box3D) {
    const extent = Vec3.sub(Vec3.zero(), box.max, box.min)

    const n = Math.pow(128, 3)
    const f = (extent[0] * extent[1] * extent[2]) / n
    const s = Math.pow(f, 1 / 3)
    const size = Vec3.zero()
    // Vec3.scale(size, extent, s)
    Vec3.ceil(size, Vec3.scale(size, extent, s))
    const delta = Vec3.div(Vec3.zero(), extent, size)
    return delta
}

async function createGaussianSurfaceMesh(ctx: RuntimeContext, unit: Unit, structure: Structure, props: GaussianSurfaceMeshProps, mesh?: Mesh): Promise<Mesh> {

    const { elements } = unit;
    const elementCount = elements.length;

    const r = 2.5;

    const v = Vec3.zero()
    const p = Vec3.zero()
    const pos = unit.conformation.invariantPosition
    const box = unit.lookup3d.boundary.box
    const expandedBox = Box3D.expand(Box3D.empty(), box, Vec3.create(r*3, r*3, r*3));
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
    const delta = getDelta(Box3D.expand(Box3D.empty(), structure.boundary.box, Vec3.create(r*3, r*3, r*3)))
    const dim = Vec3.zero()
    Vec3.ceil(dim, Vec3.mul(dim, extent, delta))
    // console.log('dim', dim, dim[0] * dim[1] * dim[2])

    const space = Tensor.Space(dim, [0, 1, 2], Float32Array)
    const data = space.create()
    const field = Tensor.create(space, data)

    for (let i = 0; i < elementCount; i++) {
        pos(elements[i], v)

        Vec3.mul(v, Vec3.sub(v, v, min), delta)

        const size = r
        const radius = size * delta[0]

        const minX = Math.floor(v[0] - radius)
        const minY = Math.floor(v[1] - radius)
        const minZ = Math.floor(v[2] - radius)
        const maxX = Math.floor(v[0] + radius)
        const maxY = Math.floor(v[1] + radius)
        const maxZ = Math.floor(v[2] + radius)

        for (let x = minX; x <= maxX; ++x) {
            for (let y = minY; y <= maxY; ++y) {
                for (let z = minZ; z <= maxZ; ++z) {
                    const dist = Vec3.distance(Vec3.set(p, x, y, z), v)
                    if (dist <= radius) {
                        const density = 1.0 - smoothstep(0.0, radius * 1.0, dist)
                        space.set(data, x, y, z, space.get(data, x, y, z) + density)
                    }
                } 
            }   
        }

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Gaussian surface', current: i, max: elementCount });
        }
    }

    // console.log('data', data)

    const surface = await computeMarchingCubes({
        isoLevel: 0.1,
        scalarField: field,
        oldSurface: mesh
        
    }).runAsChild(ctx);

    const t = Mat4.identity()
    Mat4.fromUniformScaling(t, 1 / delta[0])
    Mat4.setTranslation(t, expandedBox.min)

    ValueCell.update(surface.groupBuffer, new Float32Array(surface.vertexCount)) 
    Mesh.transformImmediate(surface, t)
    await Mesh.computeNormals(surface).runAsChild(ctx)

    // console.log('surface', surface)

    // const transform = VolumeData.getGridToCartesianTransform(volume);
    // ctx.update({ message: 'Transforming mesh...' });
    // Mesh.transformImmediate(surface, transform);

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
}
export type GaussianSurfaceProps = typeof DefaultGaussianSurfaceProps

export function GaussianSurfaceVisual(): UnitsVisual<GaussianSurfaceProps> {
    return UnitsMeshVisual<GaussianSurfaceProps>({
        defaultProps: DefaultGaussianSurfaceProps,
        createMesh: createGaussianSurfaceMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: MeshUpdateState, newProps: GaussianSurfaceProps, currentProps: GaussianSurfaceProps) => {}
    })
}