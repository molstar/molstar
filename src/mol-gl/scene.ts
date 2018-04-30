/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PointRenderable, MeshRenderable, Renderable } from './renderable'

import { ValueCell } from 'mol-util';
import { Context } from './webgl/context';

let _renderObjectId = 0;
function getNextId() {
    return _renderObjectId++ % 0x7FFFFFFF;
}

export type RenderData = { [k: string]: ValueCell<Helpers.TypedArray> }

export interface BaseRenderObject { id: number, type: string, props: {}, visible: boolean }
export interface MeshRenderObject extends BaseRenderObject { type: 'mesh', props: MeshRenderable.Props }
export interface PointRenderObject extends BaseRenderObject { type: 'point', props: PointRenderable.Props }
export type RenderObject = MeshRenderObject | PointRenderObject

export function createMeshRenderObject(props: MeshRenderable.Props): MeshRenderObject {
    return { id: getNextId(), type: 'mesh', props, visible: true }
}
export function createPointRenderObject(props: PointRenderable.Props): PointRenderObject {
    return { id: getNextId(), type: 'point', props, visible: true }
}

export function createRenderable(ctx: Context, o: RenderObject) {
    switch (o.type) {
        case 'mesh': return MeshRenderable.create(ctx, o.props)
        case 'point': return PointRenderable.create(ctx, o.props)
    }
}

interface Scene {
    add: (o: RenderObject) => void
    remove: (o: RenderObject) => void
    clear: () => void
    forEach: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => void
    count: number
}

namespace Scene {
    export function create(ctx: Context): Scene {
        const renderableMap = new Map<RenderObject, Renderable<any>>()

        return {
            add: (o: RenderObject) => {
                if (!renderableMap.has(o)) {
                    renderableMap.set(o, createRenderable(ctx, o))
                } else {
                    console.warn(`RenderObject with id '${o.id}' already present`)
                }
            },
            remove: (o: RenderObject) => {
                const renderable = renderableMap.get(o)
                if (renderable) {
                    renderable.dispose()
                    renderableMap.delete(o)
                }
            },
            clear: () => {
                renderableMap.forEach(renderable => renderable.dispose())
                renderableMap.clear()
            },
            forEach: (callbackFn: (value: Renderable<any>, key: RenderObject) => void) => {
                renderableMap.forEach(callbackFn)
            },
            get count() {
                return renderableMap.size
            }
        }
    }
}

export default Scene