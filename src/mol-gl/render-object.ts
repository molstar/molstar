/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PointRenderable, MeshRenderable, RenderableState } from './renderable'
import { RenderableValues } from './renderable/schema';
import { idFactory } from 'mol-util/id-factory';
import { Context } from './webgl/context';
import { MeshValues } from './renderable/mesh';
import { PointValues } from './renderable/point';

const getNextId = idFactory(0, 0x7FFFFFFF)

export interface BaseRenderObject { id: number, type: string, values: RenderableValues, state: RenderableState }
export interface MeshRenderObject extends BaseRenderObject { type: 'mesh', values: MeshValues }
export interface PointRenderObject extends BaseRenderObject { type: 'point', values: PointValues }
export type RenderObject = MeshRenderObject | PointRenderObject

export function createMeshRenderObject(values: MeshValues, state: RenderableState): MeshRenderObject {
    return { id: getNextId(), type: 'mesh', values, state }
}
export function createPointRenderObject(values: PointValues, state: RenderableState): PointRenderObject {
    return { id: getNextId(), type: 'point', values, state }
}

export function createRenderable(ctx: Context, o: RenderObject) {
    switch (o.type) {
        case 'mesh': return MeshRenderable(ctx, o.values, o.state)
        case 'point': return PointRenderable(ctx, o.values, o.state)
    }
}