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
import { DirectVolume2dValues, DirectVolume2dRenderable } from './renderable/direct-volume';

const getNextId = idFactory(0, 0x7FFFFFFF)

export interface BaseRenderObject { id: number, type: string, values: RenderableValues, state: RenderableState }
export interface MeshRenderObject extends BaseRenderObject { type: 'mesh', values: MeshValues }
export interface PointsRenderObject extends BaseRenderObject { type: 'points', values: PointsValues }
export interface LinesRenderObject extends BaseRenderObject { type: 'lines', values: LinesValues }
export interface GaussianDensityRenderObject extends BaseRenderObject { type: 'gaussian-density', values: GaussianDensityValues }
export interface DirectVolumeRenderObject extends BaseRenderObject { type: 'direct-volume', values: DirectVolume2dValues }
export type RenderObject = MeshRenderObject | PointsRenderObject | LinesRenderObject | GaussianDensityRenderObject | DirectVolumeRenderObject

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
export function createDirectVolumeRenderObject(values: DirectVolume2dValues, state: RenderableState): DirectVolumeRenderObject {
    return { id: getNextId(), type: 'direct-volume', values, state }
}

export function createRenderable(ctx: Context, o: RenderObject): Renderable<any> {
    switch (o.type) {
        case 'mesh': return MeshRenderable(ctx, o.id, o.values, o.state)
        case 'points': return PointsRenderable(ctx, o.id, o.values, o.state)
        case 'lines': return LinesRenderable(ctx, o.id, o.values, o.state)
        case 'gaussian-density': return GaussianDensityRenderable(ctx, o.id, o.values, o.state)
        case 'direct-volume': return DirectVolume2dRenderable(ctx, o.id, o.values, o.state)
    }
}