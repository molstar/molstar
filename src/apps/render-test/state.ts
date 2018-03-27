/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import * as glContext from 'mol-gl/context'
import { Camera } from 'mol-gl/camera'
import { Vec3 } from 'mol-math/linear-algebra'

export default class State {
    regl: REGL.Regl

    initRegl (container: HTMLDivElement) {
        const regl = glContext.create({
            container,
            extensions: [
                'OES_texture_float',
                'OES_texture_float_linear',
                // 'ext_disjoint_timer_query',
                'EXT_blend_minmax'
            ],
            // profile: true
        })

        const camera = Camera.create(regl, container, {
            center: Vec3.create(0, 0, 0)
        })

        const drawTriangle = regl({
            frag: `
            void main() {
                gl_FragColor = vec4(1, 0, 0, 1);
            }`,

            vert: `
            precision mediump float;
            uniform mat4 projection, view;
            attribute vec3 position;
            void main () {
                gl_Position = projection * view * vec4(position, 1.0);
            }`,

            attributes: {
                position: [[0, -1, 0], [-1, 0, 0], [1, 1, 0]]
            },

            count: 3
        })

        regl.frame(() => {
            camera.update((state: any) => {
                if (!camera.isDirty()) return
                console.log(state, camera)
                regl.clear({color: [0, 0, 0, 1]})
                drawTriangle()
            }, undefined)
        })

        this.regl = regl
    }

    constructor() {

    }
}
