/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PointsRenderable, MeshRenderable, RenderableState, MeshValues, PointsValues, LinesValues, LinesRenderable, Renderable } from './renderable'
import { RenderableValues } from './renderable/schema';
import { idFactory } from 'mol-util/id-factory';
import { Context } from './webgl/context';
import { GaussianDensityValues, GaussianDensityRenderable } from './renderable/gaussian-density';
import { DirectVolume2dValues, DirectVolume2dRenderable, DirectVolume3dValues, DirectVolume3dRenderable } from './renderable/direct-volume';

const getNextId = idFactory(0, 0x7FFFFFFF)

export interface BaseRenderObject { id: number, type: string, values: RenderableValues, state: RenderableState }
export interface MeshRenderObject extends BaseRenderObject { type: 'mesh', values: MeshValues }
export interface PointsRenderObject extends BaseRenderObject { type: 'points', values: PointsValues }
export interface LinesRenderObject extends BaseRenderObject { type: 'lines', values: LinesValues }
export interface GaussianDensityRenderObject extends BaseRenderObject { type: 'gaussian-density', values: GaussianDensityValues }
export interface DirectVolume2dRenderObject extends BaseRenderObject { type: 'direct-volume-2d', values: DirectVolume2dValues }
export interface DirectVolume3dRenderObject extends BaseRenderObject { type: 'direct-volume-3d', values: DirectVolume3dValues }

export type RenderObject = MeshRenderObject | PointsRenderObject | LinesRenderObject | GaussianDensityRenderObject | DirectVolume2dRenderObject | DirectVolume3dRenderObject

//

export function createMeshRenderObject(values: MeshValues, state: RenderableState): MeshRenderObject {
    return { id: getNextId(), type: 'mesh', values, state }
}
export function createPointsRenderObject(values: PointsValues, state: RenderableState): PointsRenderObject {
    return { id: getNextId(), type: 'points', values, state }
}
export function createLinesRenderObject(values: LinesValues, state: RenderableState): LinesRenderObject {
    return { id: getNextId(), type: 'lines', values, state }
}
export function createGaussianDensityRenderObject(values: GaussianDensityValues, state: RenderableState): GaussianDensityRenderObject {
    return { id: getNextId(), type: 'gaussian-density', values, state }
}
export function createDirectVolume2dRenderObject(values: DirectVolume2dValues, state: RenderableState): DirectVolume2dRenderObject {
    return { id: getNextId(), type: 'direct-volume-2d', values, state }
}
export function createDirectVolume3dRenderObject(values: DirectVolume3dValues, state: RenderableState): DirectVolume3dRenderObject {
    return { id: getNextId(), type: 'direct-volume-3d', values, state }
}

export function createRenderable(ctx: Context, o: RenderObject): Renderable<any> {
    switch (o.type) {
        case 'mesh': return MeshRenderable(ctx, o.id, o.values, o.state)
        case 'points': return PointsRenderable(ctx, o.id, o.values, o.state)
        case 'lines': return LinesRenderable(ctx, o.id, o.values, o.state)
        case 'gaussian-density': return GaussianDensityRenderable(ctx, o.id, o.values, o.state)
        case 'direct-volume-2d': return DirectVolume2dRenderable(ctx, o.id, o.values, o.state)
        case 'direct-volume-3d': return DirectVolume3dRenderable(ctx, o.id, o.values, o.state)
    }
}