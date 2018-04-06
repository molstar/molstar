/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import * as glContext from './context'
import { Camera } from './camera'
import { PointRenderable, MeshRenderable, Renderable } from './renderable'

import { Vec3, Mat4 } from 'mol-math/linear-algebra'

let _renderObjectId = 0;
function getNextId() {
    return _renderObjectId++ % 0x7FFFFFFF;
}

export interface RenderObject {
    id: number
    type: 'mesh' | 'point'
    data: any
    elements: any
    uniforms: any
}

export function createRenderObject(type: 'mesh' | 'point', data: any, elements?: any, uniforms?: any) {
    return { id: getNextId(), type, data, elements, uniforms }
}

export interface Renderer {
    add: (o: RenderObject) => void
    remove: (o: RenderObject) => void
    draw: () => void
    frame: () => void
}

export function createRenderable(regl: REGL.Regl, o: RenderObject) {
    switch (o.type) {
        case 'mesh': return MeshRenderable.create(regl, o.data, o.uniforms || {}, o.elements)
        case 'point': return PointRenderable.create(regl, o.data)
    }
}

export function createRenderer(container: HTMLDivElement): Renderer {
    const renderableList: Renderable<any>[] = []
    const objectIdRenderableMap: { [k: number]: Renderable<any> } = {}

    const regl = glContext.create({
        container,
        extensions: [
            'OES_texture_float',
            'OES_texture_float_linear',
            'OES_element_index_uint',
            // 'EXT_disjoint_timer_query',
            'EXT_blend_minmax',
            'ANGLE_instanced_arrays'
        ],
        // profile: true
    })

    const camera = Camera.create(regl, container, {
        center: Vec3.create(0, 0, 0),
        near: 0.01,
        far: 1000
    })

    const baseContext = regl({
        context: {
            model: Mat4.identity(),
            transform: Mat4.setTranslation(Mat4.identity(), Vec3.create(6, 0, 0))
        },
        uniforms: {
            model: regl.context('model' as any),
            transform: regl.context('transform' as any),
            'light.position': Vec3.create(0, 0, -100),
            'light.color': Vec3.create(1.0, 1.0, 1.0),
            'light.ambient': Vec3.create(0.5, 0.5, 0.5),
            'light.falloff': 0,
            'light.radius': 500
        }
    })

    const draw = () => {
        camera.update((state: any) => {
            if (!camera.isDirty()) return;
            baseContext(() => {
                // console.log(ctx)
                regl.clear({color: [0, 0, 0, 1]})
                // TODO painters sort, filter visible, filter picking, visibility culling?
                renderableList.forEach(r => {
                    r.draw()
                })
            })
        }, undefined)
    }

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
        draw,
        frame: () => {
            regl.frame((ctx) => draw())
        }
    }
}