/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createGl } from './gl.shim';

import { PerspectiveCamera } from 'mol-view/camera/perspective';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { ValueCell } from 'mol-util';

import Renderer from '../renderer';
import { fillSerial } from '../renderable/util';
import { createUniformColor } from 'mol-geo/util/color-data';
import { createUniformSize } from 'mol-geo/util/size-data';
import { createContext } from '../webgl/context';
import { RenderableState } from '../renderable';
import { createPointRenderObject } from '../render-object';
import { PointValues } from '../renderable/point';
import Scene from '../scene';

// function writeImage(gl: WebGLRenderingContext, width: number, height: number) {
//     const pixels = new Uint8Array(width * height * 4)
//     gl.readPixels(0, 0, width, height, gl.RGBA, gl.UNSIGNED_BYTE, pixels)
//     process.stdout.write(['P3\n# gl.ppm\n', width, ' ', height, '\n255\n'].join(''))
//     for (let i = 0; i<pixels.length; i+=4) {
//         for (let j = 0; j<3; ++j) {
//             process.stdout.write(pixels[i+j] + ' ')
//         }
//     }
// }

function createRenderer(gl: WebGLRenderingContext) {
    const ctx = createContext(gl)
    const camera = PerspectiveCamera.create({
        near: 0.01,
        far: 10000,
        position: Vec3.create(0, 0, 50)
    })
    const renderer = Renderer.create(ctx, camera)
    return { ctx, camera, renderer }
}

function createPoints() {
    const aPosition = ValueCell.create(new Float32Array([0, -1, 0, -1, 0, 0, 1, 1, 0]))
    const aElementId = ValueCell.create(fillSerial(new Float32Array(3)))
    const aInstanceId = ValueCell.create(fillSerial(new Float32Array(1)))
    const color = createUniformColor({ value: 0xFF0000 })
    const size = createUniformSize({ value: 1 })

    const aTransform = ValueCell.create(new Float32Array(16))
    const m4 = Mat4.identity()
    Mat4.toArray(m4, aTransform.ref.value, 0)

    const values: PointValues = {
        aPosition,
        aElementId,
        aTransform,
        aInstanceId,
        ...color,
        ...size,

        uAlpha: ValueCell.create(1.0),
        uInstanceCount: ValueCell.create(1),
        uElementCount: ValueCell.create(3),

        drawCount: ValueCell.create(3),
        instanceCount: ValueCell.create(1),

        dPointSizeAttenuation: ValueCell.create(true)
    }
    const state: RenderableState = {
        visible: true,
        depthMask: true,
    }

    return createPointRenderObject(values, state)
}

describe('renderer', () => {
    it('basic', () => {
        const [ width, height ] = [ 32, 32 ]
        const gl = createGl(width, height, { preserveDrawingBuffer: true })
        const { ctx, renderer } = createRenderer(gl)

        expect(ctx.gl.canvas.width).toBe(32)
        expect(ctx.gl.canvas.height).toBe(32)

        expect(ctx.bufferCount).toBe(0);
        expect(ctx.textureCount).toBe(0);
        expect(ctx.vaoCount).toBe(0);
        expect(ctx.programCache.count).toBe(0);
        expect(ctx.shaderCache.count).toBe(0);

        renderer.setViewport({ x: 0, y: 0, width: 64, height: 48 })
        expect(ctx.gl.getParameter(ctx.gl.VIEWPORT)[2]).toBe(64)
        expect(ctx.gl.getParameter(ctx.gl.VIEWPORT)[3]).toBe(48)
    })

    it('points', () => {
        const [ width, height ] = [ 32, 32 ]
        const gl = createGl(width, height, { preserveDrawingBuffer: true })
        const { ctx } = createRenderer(gl)
        const scene = Scene.create(ctx)

        const points = createPoints()

        scene.add(points)
        expect(ctx.bufferCount).toBe(6);
        expect(ctx.textureCount).toBe(1);
        expect(ctx.vaoCount).toBe(1);
        expect(ctx.programCache.count).toBe(1);
        expect(ctx.shaderCache.count).toBe(2);

        scene.remove(points)
        expect(ctx.bufferCount).toBe(0);
        expect(ctx.textureCount).toBe(0);
        expect(ctx.vaoCount).toBe(0);
        expect(ctx.programCache.count).toBe(1);
        expect(ctx.shaderCache.count).toBe(2);

        ctx.programCache.dispose()
        expect(ctx.programCache.count).toBe(0);

        ctx.shaderCache.clear()
        expect(ctx.shaderCache.count).toBe(0);
        // console.log('moin', ctx)
    })
})