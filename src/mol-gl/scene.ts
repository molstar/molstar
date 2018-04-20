/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { PointRenderable, MeshRenderable, Renderable } from './renderable'

import { ValueCell } from 'mol-util';

let _renderObjectId = 0;
function getNextId() {
    return _renderObjectId++ % 0x7FFFFFFF;
}

export type RenderData = { [k: string]: ValueCell<Helpers.TypedArray> }

export interface MeshRenderObject { id: number, type: 'mesh', props: MeshRenderable.Data }
export interface PointRenderObject { id: number, type: 'point', props: PointRenderable.Data }
export type RenderObject = MeshRenderObject | PointRenderObject

export function createMeshRenderObject(props: MeshRenderable.Data): MeshRenderObject {
    return { id: getNextId(), type: 'mesh', props }
}
export function createPointRenderObject(props: PointRenderable.Data): PointRenderObject {
    return { id: getNextId(), type: 'point', props }
}

export function createRenderable(regl: REGL.Regl, o: RenderObject) {
    switch (o.type) {
        case 'mesh': return MeshRenderable.create(regl, o.props)
        case 'point': return PointRenderable.create(regl, o.props)
    }
}

interface Scene {
    add: (o: RenderObject) => void
    remove: (o: RenderObject) => void
    clear: () => void
    forEach: (callbackFn: (value: Renderable) => void) => void
    count: number
}

namespace Scene {
    export function create(regl: REGL.Regl): Scene {
        const renderableList: Renderable[] = []
        const objectIdRenderableMap: { [k: number]: Renderable } = {}

        return {
            add: (o: RenderObject) => {
                const renderable = createRenderable(regl, o)
                renderableList.push(renderable)
                objectIdRenderableMap[o.id] = renderable
            },
            remove: (o: RenderObject) => {
                if (o.id in objectIdRenderableMap) {
                    objectIdRenderableMap[o.id].dispose()
                    delete objectIdRenderableMap[o.id]
                }
            },
            clear: () => {
                for (const id in objectIdRenderableMap) {
                    objectIdRenderableMap[id].dispose()
                    delete objectIdRenderableMap[id]
                }
                renderableList.length = 0
            },
            forEach: (callbackFn: (value: Renderable) => void) => {
                renderableList.forEach(callbackFn)
            },
            get count() {
                return renderableList.length
            }
        }
    }
}

export default Scene