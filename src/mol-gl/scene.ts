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

export interface RenderObject {
    id: number
    type: 'mesh' | 'point'
    data: PointRenderable.Data | MeshRenderable.Data
    uniforms: { [k: string]: REGL.Uniform }
}

export function createRenderObject(type: 'mesh' | 'point', data: PointRenderable.Data | MeshRenderable.Data, uniforms: { [k: string]: REGL.Uniform }) {
    return { id: getNextId(), type, data, uniforms }
}

export function createRenderable(regl: REGL.Regl, o: RenderObject) {
    switch (o.type) {
        case 'mesh': return MeshRenderable.create(regl, o.data as MeshRenderable.Data, o.uniforms || {})
        case 'point': return PointRenderable.create(regl, o.data as PointRenderable.Data)
    }
}

interface Scene {
    add: (o: RenderObject) => void
    remove: (o: RenderObject) => void
    clear: () => void
    forEach: (callbackFn: (value: Renderable) => void) => void
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
                    // TODO
                    // objectIdRenderableMap[o.id].destroy()
                    delete objectIdRenderableMap[o.id]
                }
            },
            clear: () => {
                for (const id in objectIdRenderableMap) {
                    // TODO
                    // objectIdRenderableMap[id].destroy()
                    delete objectIdRenderableMap[id]
                }
                renderableList.length = 0
            },
            forEach: (callbackFn: (value: Renderable) => void) => {
                renderableList.forEach(callbackFn)
            }
        }
    }
}

export default Scene