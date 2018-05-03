/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData, VolumeIsoValue } from 'mol-model/volume'
import { Task } from 'mol-task'
import { computeMarchingCubes } from '../../util/marching-cubes/algorithm';
import { Mesh } from '../../shape/mesh';
import { VolumeElementRepresentation } from '.';
import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/scene';
import { fillSerial } from 'mol-gl/renderable/util';
import { ValueCell } from 'mol-util';
import { Mat4 } from 'mol-math/linear-algebra';
import { createUniformColor } from '../../util/color-data';

export function computeVolumeSurface(volume: VolumeData, isoValue: VolumeIsoValue) {
    return Task.create<Mesh>('Volume Surface', async ctx => {
        ctx.update({ message: 'Marching cubes...' });

        const mesh = await ctx.runChild(computeMarchingCubes({
            isoLevel: VolumeIsoValue.toAbsolute(isoValue).absoluteValue,
            scalarField: volume.data
        }));

        const transform = VolumeData.getGridToCartesianTransform(volume);
        ctx.update({ message: 'Transforming mesh...' });
        Mesh.transformImmediate(mesh, transform);

        return mesh;
    });
}

export const DefaultSurfaceProps = {
    isoValue: VolumeIsoValue.relative({ min: 0, max: 0, mean: 0, sigma: 0 }, 0),
    alpha: 0.5,
    flatShaded: true,
    flipSided: true,
    doubleSided: true
}
export type SurfaceProps = Partial<typeof DefaultSurfaceProps>

export default function Surface(): VolumeElementRepresentation<SurfaceProps> {
    const renderObjects: RenderObject[] = []
    let surface: MeshRenderObject
    let curProps = DefaultSurfaceProps

    return {
        renderObjects,
        create(volume: VolumeData, props: SurfaceProps = {}) {
            return Task.create('Point.create', async ctx => {
                renderObjects.length = 0 // clear
                curProps = { ...DefaultSurfaceProps, ...props }
                const { alpha, flatShaded, flipSided, doubleSided } = curProps

                const mesh = await ctx.runChild(computeVolumeSurface(volume, curProps.isoValue))

                surface = createMeshRenderObject({
                    objectId: 0,
                    alpha,

                    position: mesh.vertexBuffer,
                    normal: mesh.normalBuffer,
                    id: ValueCell.create(fillSerial(new Float32Array(mesh.vertexCount / 3))),
                    color: createUniformColor({ value: 0x7ec0ee }),
                    transform: ValueCell.create(new Float32Array(Mat4.identity())),
                    index: mesh.indexBuffer,

                    instanceCount: 1,
                    indexCount: mesh.triangleCount,
                    elementCount: mesh.triangleCount,
                    positionCount: mesh.vertexCount / 3,

                    flatShaded,
                    doubleSided,
                    flipSided
                })
                renderObjects.push(surface)
            })
        },
        update(props: SurfaceProps) {
            return Task.create('Surface.update', async ctx => {
                // TODO
                return false
            })
        }
    }
}
