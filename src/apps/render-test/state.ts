/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import * as glContext from 'mol-gl/context'
import { Camera } from 'mol-gl/camera'
import { Vec3 } from 'mol-math/linear-algebra'
import Point from 'mol-gl/renderable/point'

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

        const points = Point.create(regl, {
            position: new Float32Array([0, -1, 0, -1, 0, 0, 1, 1, 0])
        })

        regl.frame(() => {
            camera.update((state: any) => {
                if (!camera.isDirty()) return
                regl.clear({color: [0, 0, 0, 1]})
                points.update(a => { a.position[0] = Math.random() })
                points.draw()
            }, undefined)
        })

        this.regl = regl
    }

    constructor() {

    }
}
