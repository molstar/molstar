/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import createGl = require('gl');

import { PerspectiveCamera } from 'mol-view/camera/perspective';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { ValueCell } from 'mol-util';

import Renderer from '../renderer';
import { createPointRenderObject } from '../scene';
import { fillSerial } from '../renderable/util';
import { createUniformColor } from 'mol-geo/util/color-data';
import { createUniformSize } from 'mol-geo/util/size-data';
import { createContext } from '../webgl/context';

function writeImage(gl: WebGLRenderingContext, width: number, height: number) {
    const pixels = new Uint8Array(width * height * 4)
    gl.readPixels(0, 0, width, height, gl.RGBA, gl.UNSIGNED_BYTE, pixels)
    process.stdout.write(['P3\n# gl.ppm\n', width, ' ', height, '\n255\n'].join(''))
    for (let i = 0; i<pixels.length; i+=4) {
        for (let j = 0; j<3; ++j) {
            process.stdout.write(pixels[i+j] + ' ')
        }
    }
}

function createRenderer(gl: WebGLRenderingContext) {
    const ctx = createContext(gl)
    const camera = PerspectiveCamera.create({
        near: 0.01,
        far: 10000,
        position: Vec3.create(0, 0, 50)
    })
    return Renderer.create(ctx, camera)
}

function createPoints() {
    const position = ValueCell.create(new Float32Array([0, -1, 0, -1, 0, 0, 1, 1, 0]))
    const id = ValueCell.create(fillSerial(new Float32Array(3)))
    const color = createUniformColor({ value: 0xFF0000 })
    const size = createUniformSize({ value: 1 })

    const transform = ValueCell.create(new Float32Array(16))
    const m4 = Mat4.identity()
    Mat4.toArray(m4, transform.ref.value, 0)

    return createPointRenderObject({
        objectId: 0,
        alpha: 1.0,
        visible: true,

        position,
        id,
        color,
        size,
        transform,

        instanceCount: 1,
        elementCount: 3,
        positionCount: 3
    })
}

// TODO not working
// - shaders not transformed via glslify
describe.skip('renderer', () => {
    it('basic', () => {
        const [ width, height ] = [ 32, 32 ]
        const gl = createGl(width, height, { preserveDrawingBuffer: true })
        const renderer = createRenderer(gl)
        const points = createPoints()

        renderer.add(points)
        renderer.draw()

        writeImage(gl, width, height)
    })
})