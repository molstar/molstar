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

const getNextId = idFactory(0, 0x7FFFFFFF)

export const getNextMaterialId = idFactory(0, 0x7FFFFFFF)

export interface BaseRenderObject { id: number, type: string, values: RenderableValues, state: RenderableState, materialId: number }
export interface MeshRenderObject extends BaseRenderObject { type: 'mesh', values: MeshValues }
export interface PointsRenderObject extends BaseRenderObject { type: 'points', values: PointsValues }
export interface SpheresRenderObject extends BaseRenderObject { type: 'spheres', values: SpheresValues }
export interface TextRenderObject extends BaseRenderObject { type: 'text', values: TextValues }
export interface LinesRenderObject extends BaseRenderObject { type: 'lines', values: LinesValues }
export interface DirectVolumeRenderObject extends BaseRenderObject { type: 'direct-volume', values: DirectVolumeValues }

//

export type GraphicsRenderObject = MeshRenderObject | PointsRenderObject | SpheresRenderObject | TextRenderObject | LinesRenderObject | DirectVolumeRenderObject

export type RenderObjectKindType = {
    'mesh': MeshRenderObject
    'points': PointsRenderObject
    'spheres': SpheresRenderObject
    'text': TextRenderObject
    'lines': LinesRenderObject
    'direct-volume': DirectVolumeRenderObject
}
export type RenderObjectValuesType = {
    'mesh': MeshValues
    'points': PointsValues
    'spheres': SpheresValues
    'text': TextValues
    'lines': LinesValues
    'direct-volume': DirectVolumeValues
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
    }
}