/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import * as glContext from 'mol-gl/context'
import { Camera } from 'mol-gl/camera'
import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import PointRenderable from 'mol-gl/renderable/point'
import MeshRenderable from 'mol-gl/renderable/mesh'
import Attribute from 'mol-gl/attribute';
import Model from 'mol-gl/model';
// import { positionFromModel } from 'mol-geo/shape/point'

export default class State {
    regl: REGL.Regl

    initRegl (container: HTMLDivElement) {
        const regl = glContext.create({
            container,
            extensions: [
                'OES_texture_float',
                'OES_texture_float_linear',
                'OES_element_index_uint',
                // 'ext_disjoint_timer_query',
                'EXT_blend_minmax',
                'ANGLE_instanced_arrays'
            ],
            // profile: true
        })

        const camera = Camera.create(regl, container, {
            center: Vec3.create(0, 0, 0)
        })

        const p1 = Vec3.create(0, 4, 0)
        const p2 = Vec3.create(-3, 0, 0)

        const model1 = Model(regl)
        const model2 = Model(regl, { position: p1 })
        const model3 = Model(regl, { position: p2 })

        const position = Attribute.create(regl, new Float32Array([0, -1, 0, -1, 0, 0, 1, 1, 0]), 3)

        const transformArray1 = new Float32Array(16)
        const transformArray2 = new Float32Array(16 * 3)
        const m4 = Mat4.identity()
        Mat4.toArray(m4, transformArray1)
        Mat4.toArray(m4, transformArray2)
        Mat4.setTranslation(m4, p1)
        Mat4.toArray(m4, transformArray2, 16)
        Mat4.setTranslation(m4, p2)
        Mat4.toArray(m4, transformArray2, 32)
        const transform1 = Attribute.create(regl, transformArray1, 16, 1)
        const transform2 = Attribute.create(regl, transformArray2, 16, 1)

        // TODO use https://github.com/substack/glsl-matrix-texture

        // position.update((array: Float32Array) => {
        //     positionFromModel({}, array, 0)
        // })

        const points = PointRenderable.create(regl, { position, transform: transform1 })
        const mesh = MeshRenderable.create(regl, { position, transform: transform2 })

        const baseContext = regl({
            context: {
                model: Mat4.identity(),
                transform: Mat4.setTranslation(Mat4.identity(), Vec3.create(6, 0, 0))
            },
            uniforms: {
                model: regl.context('model' as any),
                transform: regl.context('transform' as any),
            }
        })

        regl.frame((ctx) => {
            camera.update((state: any) => {
                if (!camera.isDirty()) return
                baseContext(() => {
                    console.log(ctx)
                    regl.clear({color: [0, 0, 0, 1]})
                    position.update(array => { array[0] = Math.random() })
                    // points.update(a => { a.position[0] = Math.random() })
                    mesh.draw()
                    model1({}, ({ transform }) => {
                        points.draw()
                    })
                    model2({}, ({ transform }) => {
                        points.draw()
                        model3({ transform }, () => {
                            points.draw()
                        })
                    })
                })
            }, undefined)
        })

        this.regl = regl
    }

    constructor() {

    }
}
