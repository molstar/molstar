/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RenderableState, Renderable } from './renderable'
import { RenderableValues } from './renderable/schema';
import { idFactory } from 'mol-util/id-factory';
import { WebGLContext } from './webgl/context';
import { DirectVolumeValues, DirectVolumeRenderable } from './renderable/direct-volume';
import { MeshValues, MeshRenderable } from './renderable/mesh';
import { PointsValues, PointsRenderable } from './renderable/points';
import { LinesValues, LinesRenderable } from './renderable/lines';
import { SpheresValues, SpheresRenderable } from './renderable/spheres';
import { TextValues, TextRenderable } from './renderable/text';
import { TextureMeshValues, TextureMeshRenderable } from './renderable/texture-mesh';

const getNextId = idFactory(0, 0x7FFFFFFF)

export const getNextMaterialId = idFactory(0, 0x7FFFFFFF)

export interface BaseRenderObject<T extends RenderableValues> { id: number, type: string, values: T, state: RenderableState, materialId: number }
export interface MeshRenderObject extends BaseRenderObject<MeshValues> { type: 'mesh' }
export interface PointsRenderObject extends BaseRenderObject<PointsValues> { type: 'points' }
export interface SpheresRenderObject extends BaseRenderObject<SpheresValues> { type: 'spheres' }
export interface TextRenderObject extends BaseRenderObject<TextValues> { type: 'text' }
export interface LinesRenderObject extends BaseRenderObject<LinesValues> { type: 'lines' }
export interface DirectVolumeRenderObject extends BaseRenderObject<DirectVolumeValues> { type: 'direct-volume' }
export interface TextureMeshRenderObject extends BaseRenderObject<TextureMeshValues> { type: 'texture-mesh' }

//

export type GraphicsRenderObject = MeshRenderObject | PointsRenderObject | SpheresRenderObject | TextRenderObject | LinesRenderObject | DirectVolumeRenderObject | TextureMeshRenderObject

export type RenderObjectKindType = {
    'mesh': MeshRenderObject
    'points': PointsRenderObject
    'spheres': SpheresRenderObject
    'text': TextRenderObject
    'lines': LinesRenderObject
    'direct-volume': DirectVolumeRenderObject
    'texture-mesh': TextureMeshRenderObject
}
export type RenderObjectValuesType = {
    'mesh': MeshValues
    'points': PointsValues
    'spheres': SpheresValues
    'text': TextValues
    'lines': LinesValues
    'direct-volume': DirectVolumeValues
    'texture-mesh': TextureMeshValues
}
export type RenderObjectType = keyof RenderObjectKindType

//

export function createRenderObject<T extends RenderObjectType>(type: T, values: RenderObjectValuesType[T], state: RenderableState, materialId: number): RenderObjectKindType[T] {
    return { id: getNextId(), type, values, state, materialId } as RenderObjectKindType[T]
}

export function createRenderable(ctx: WebGLContext, o: GraphicsRenderObject): Renderable<any> {
    switch (o.type) {
        case 'mesh': return MeshRenderable(ctx, o.id, o.values, o.state, o.materialId)
        case 'points': return PointsRenderable(ctx, o.id, o.values, o.state, o.materialId)
        case 'spheres': return SpheresRenderable(ctx, o.id, o.values, o.state, o.materialId)
        case 'text': return TextRenderable(ctx, o.id, o.values, o.state, o.materialId)
        case 'lines': return LinesRenderable(ctx, o.id, o.values, o.state, o.materialId)
        case 'direct-volume': return DirectVolumeRenderable(ctx, o.id, o.values, o.state, o.materialId)
        case 'texture-mesh': return TextureMeshRenderable(ctx, o.id, o.values, o.state, o.materialId)
    }
}