/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RenderableState, Renderable } from './renderable';
import { idFactory } from '../mol-util/id-factory';
import { WebGLContext } from './webgl/context';
import { DirectVolumeValues, DirectVolumeRenderable } from './renderable/direct-volume';
import { MeshValues, MeshRenderable } from './renderable/mesh';
import { PointsValues, PointsRenderable } from './renderable/points';
import { LinesValues, LinesRenderable } from './renderable/lines';
import { SpheresValues, SpheresRenderable } from './renderable/spheres';
import { TextValues, TextRenderable } from './renderable/text';
import { TextureMeshValues, TextureMeshRenderable } from './renderable/texture-mesh';
import { ImageValues, ImageRenderable } from './renderable/image';

const getNextId = idFactory(0, 0x7FFFFFFF);

export const getNextMaterialId = idFactory(0, 0x7FFFFFFF);

export interface GraphicsRenderObject<T extends RenderObjectType = RenderObjectType> {
    readonly id: number,
    readonly type: T,
    readonly values: RenderObjectValues<T>,
    readonly state: RenderableState,
    readonly materialId: number
}

export type RenderObjectType = 'mesh' | 'points' | 'spheres' | 'text' | 'lines' | 'direct-volume' | 'image' | 'texture-mesh'

export type RenderObjectValues<T extends RenderObjectType> =
    T extends 'mesh' ? MeshValues :
        T extends 'points' ? PointsValues :
            T extends 'spheres' ? SpheresValues :
                T extends 'text' ? TextValues :
                    T extends 'lines' ? LinesValues :
                        T extends 'direct-volume' ? DirectVolumeValues :
                            T extends 'image' ? ImageValues :
                                T extends 'texture-mesh' ? TextureMeshValues : never

//

export function createRenderObject<T extends RenderObjectType>(type: T, values: RenderObjectValues<T>, state: RenderableState, materialId: number): GraphicsRenderObject<T> {
    return { id: getNextId(), type, values, state, materialId } as GraphicsRenderObject<T>;
}

export function createRenderable<T extends RenderObjectType>(ctx: WebGLContext, o: GraphicsRenderObject<T>): Renderable<any> {
    switch (o.type) {
        case 'mesh': return MeshRenderable(ctx, o.id, o.values as MeshValues, o.state, o.materialId);
        case 'points': return PointsRenderable(ctx, o.id, o.values as PointsValues, o.state, o.materialId);
        case 'spheres': return SpheresRenderable(ctx, o.id, o.values as SpheresValues, o.state, o.materialId);
        case 'text': return TextRenderable(ctx, o.id, o.values as TextValues, o.state, o.materialId);
        case 'lines': return LinesRenderable(ctx, o.id, o.values as LinesValues, o.state, o.materialId);
        case 'direct-volume': return DirectVolumeRenderable(ctx, o.id, o.values as DirectVolumeValues, o.state, o.materialId);
        case 'image': return ImageRenderable(ctx, o.id, o.values as ImageValues, o.state, o.materialId);
        case 'texture-mesh': return TextureMeshRenderable(ctx, o.id, o.values as TextureMeshValues, o.state, o.materialId);
    }
    throw new Error('unsupported type');
}
