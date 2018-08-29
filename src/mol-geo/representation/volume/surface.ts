/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData, VolumeIsoValue } from 'mol-model/volume'
import { Task, RuntimeContext } from 'mol-task'
import { computeMarchingCubes } from '../../util/marching-cubes/algorithm';
import { Mesh } from '../../mesh/mesh';
import { VolumeVisual } from '.';
import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object';
import { ValueCell, defaults } from 'mol-util';
import { Mat4 } from 'mol-math/linear-algebra';
import { createValueColor } from '../../util/color-data';
import { getMeshData } from '../../util/mesh-data';
import { RenderableState, MeshValues } from 'mol-gl/renderable';
import { PickingId } from '../../util/picking';
import { createEmptyMarkers, MarkerAction } from '../../util/marker-data';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { fillSerial } from 'mol-util/array';
import { Color } from 'mol-util/color';

export function computeVolumeSurface(volume: VolumeData, isoValue: VolumeIsoValue) {
    return Task.create<Mesh>('Volume Surface', async ctx => {
        ctx.update({ message: 'Marching cubes...' });

        const mesh = await computeMarchingCubes({
            isoLevel: VolumeIsoValue.toAbsolute(isoValue).absoluteValue,
            scalarField: volume.data
        }).runAsChild(ctx);

        const transform = VolumeData.getGridToCartesianTransform(volume);
        ctx.update({ message: 'Transforming mesh...' });
        Mesh.transformImmediate(mesh, transform);

        return mesh;
    });
}

export const DefaultSurfaceProps = {
    isoValue: VolumeIsoValue.relative({ min: 0, max: 0, mean: 0, sigma: 0 }, 0),
    alpha: 0.5,
    visible: true,
    flatShaded: true,
    flipSided: true,
    doubleSided: true,
    depthMask: true,
    useFog: true
}
export type SurfaceProps = Partial<typeof DefaultSurfaceProps>

export default function SurfaceVisual(): VolumeVisual<SurfaceProps> {
    let renderObject: MeshRenderObject
    let curProps = DefaultSurfaceProps

    return {
        get renderObject () { return renderObject },
        async createOrUpdate(ctx: RuntimeContext, props: SurfaceProps = {}, volume?: VolumeData) {
            props = { ...DefaultSurfaceProps, ...props }

            if (!volume) return

            const mesh = await computeVolumeSurface(volume, curProps.isoValue).runAsChild(ctx)
            if (!props.flatShaded) {
                Mesh.computeNormalsImmediate(mesh)
            }

            const instanceCount = 1
            const color = createValueColor(Color(0x7ec0ee))
            const marker = createEmptyMarkers()

            const values: MeshValues = {
                ...getMeshData(mesh),
                aTransform: ValueCell.create(new Float32Array(Mat4.identity())),
                aInstance: ValueCell.create(fillSerial(new Float32Array(instanceCount))),
                ...color,
                ...marker,

                uAlpha: ValueCell.create(defaults(props.alpha, 1.0)),
                uInstanceCount: ValueCell.create(instanceCount),
                uGroupCount: ValueCell.create(mesh.triangleCount),

                elements: mesh.indexBuffer,

                drawCount: ValueCell.create(mesh.triangleCount * 3),
                instanceCount: ValueCell.create(instanceCount),

                dDoubleSided: ValueCell.create(defaults(props.doubleSided, true)),
                dFlatShaded: ValueCell.create(defaults(props.flatShaded, true)),
                dFlipSided: ValueCell.create(false),
                dUseFog: ValueCell.create(defaults(props.useFog, true)),
            }
            const state: RenderableState = {
                depthMask: defaults(props.depthMask, true),
                visible: defaults(props.visible, true)
            }

            renderObject = createMeshRenderObject(values, state)
        },
        getLoci(pickingId: PickingId) {
            // TODO
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            // TODO
        },
        destroy() {
            // TODO
        }
    }
}
