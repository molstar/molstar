/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RenderableState, Renderable } from './renderable'
import { RenderableValues } from './renderable/schema';
import { idFactory } from 'mol-util/id-factory';
import { WebGLContext } from './webgl/context';
import { GaussianDensityValues, GaussianDensityRenderable } from './renderable/gaussian-density';
import { DirectVolumeValues, DirectVolumeRenderable } from './renderable/direct-volume';
import { MeshValues, MeshRenderable } from './renderable/mesh';
import { PointsValues, PointsRenderable } from './renderable/points';
import { LinesValues, LinesRenderable } from './renderable/lines';
import { SpheresValues, SpheresRenderable } from './renderable/spheres';
import { TextValues, TextRenderable } from './renderable/text';

const getNextId = idFactory(0, 0x7FFFFFFF)

export interface BaseRenderObject { id: number, type: string, values: RenderableValues, state: RenderableState }
export interface MeshRenderObject extends BaseRenderObject { type: 'mesh', values: MeshValues }
export interface PointsRenderObject extends BaseRenderObject { type: 'points', values: PointsValues }
export interface SpheresRenderObject extends BaseRenderObject { type: 'spheres', values: SpheresValues }
export interface TextRenderObject extends BaseRenderObject { type: 'text', values: TextValues }
export interface LinesRenderObject extends BaseRenderObject { type: 'lines', values: LinesValues }
export interface DirectVolumeRenderObject extends BaseRenderObject { type: 'direct-volume', values: DirectVolumeValues }

export interface GaussianDensityRenderObject extends BaseRenderObject { type: 'gaussian-density', values: GaussianDensityValues }

//

export type GraphicsRenderObject = MeshRenderObject | PointsRenderObject | SpheresRenderObject | TextRenderObject | LinesRenderObject | DirectVolumeRenderObject

export type ComputeRenderObject = GaussianDensityRenderObject

export type RenderObject = GraphicsRenderObject | ComputeRenderObject

export type RenderObjectKindType = {
    'mesh': MeshRenderObject
    'points': PointsRenderObject
    'spheres': SpheresRenderObject
    'text': TextRenderObject
    'lines': LinesRenderObject
    'direct-volume': DirectVolumeRenderObject

    'gaussian-density': GaussianDensityRenderObject
}
export type RenderObjectValuesType = {
    'mesh': MeshValues
    'points': PointsValues
    'spheres': SpheresValues
    'text': TextValues
    'lines': LinesValues
    'direct-volume': DirectVolumeValues

    'gaussian-density': GaussianDensityValues
}
export type RenderObjectType = keyof RenderObjectKindType

//

export function createRenderObject<T extends RenderObjectType>(type: T, values: RenderObjectValuesType[T], state: RenderableState): RenderObjectKindType[T] {
    return { id: getNextId(), type, values, state } as RenderObjectKindType[T]
}

export function createRenderable(ctx: WebGLContext, o: RenderObject): Renderable<any> {
    switch (o.type) {
        case 'mesh': return MeshRenderable(ctx, o.id, o.values, o.state)
        case 'points': return PointsRenderable(ctx, o.id, o.values, o.state)
        case 'spheres': return SpheresRenderable(ctx, o.id, o.values, o.state)
        case 'text': return TextRenderable(ctx, o.id, o.values, o.state)
        case 'lines': return LinesRenderable(ctx, o.id, o.values, o.state)
        case 'direct-volume': return DirectVolumeRenderable(ctx, o.id, o.values, o.state)

        case 'gaussian-density': return GaussianDensityRenderable(ctx, o.id, o.values, o.state)
    }
}