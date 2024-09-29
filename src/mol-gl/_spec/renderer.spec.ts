/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createGl } from './gl.shim';
import { Camera } from '../../mol-canvas3d/camera';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Renderer } from '../renderer';
import { createContext } from '../webgl/context';
import { Scene } from '../scene';
import { createPoints } from './points.spec';

export function createRenderer(gl: WebGLRenderingContext) {
    const ctx = createContext(gl);
    const camera = new Camera({
        position: Vec3.create(0, 0, 50)
    });
    const renderer = Renderer.create(ctx);
    return { ctx, camera, renderer };
}

describe('renderer', () => {
    it('basic', () => {
        const [width, height] = [32, 32];
        const gl = createGl(width, height, { preserveDrawingBuffer: true });
        const { ctx, renderer } = createRenderer(gl);

        expect(ctx.gl.drawingBufferWidth).toBe(32);
        expect(ctx.gl.drawingBufferHeight).toBe(32);

        expect(ctx.stats.resourceCounts.attribute).toBe(0);
        expect(ctx.stats.resourceCounts.texture).toBe(1);
        expect(ctx.stats.resourceCounts.vertexArray).toBe(0);
        expect(ctx.stats.resourceCounts.program).toBe(0);
        expect(ctx.stats.resourceCounts.shader).toBe(0);

        renderer.setViewport(0, 0, 64, 48);
        expect(ctx.gl.getParameter(ctx.gl.VIEWPORT)[2]).toBe(64);
        expect(ctx.gl.getParameter(ctx.gl.VIEWPORT)[3]).toBe(48);
    });

    it('points', async () => {
        const [width, height] = [32, 32];
        const gl = createGl(width, height, { preserveDrawingBuffer: true });
        const { ctx } = createRenderer(gl);
        const scene = Scene.create(ctx);

        const points = createPoints();

        scene.add(points);
        scene.commit();
        expect(ctx.stats.resourceCounts.attribute).toBe(ctx.isWebGL2 ? 4 : 5);
        expect(ctx.stats.resourceCounts.texture).toBe(10);
        expect(ctx.stats.resourceCounts.vertexArray).toBe(ctx.extensions.vertexArrayObject ? 5 : 0);
        expect(ctx.stats.resourceCounts.program).toBe(5);
        expect(ctx.stats.resourceCounts.shader).toBe(10);

        scene.remove(points);
        scene.commit();
        expect(ctx.stats.resourceCounts.attribute).toBe(0);
        expect(ctx.stats.resourceCounts.texture).toBe(1);
        expect(ctx.stats.resourceCounts.vertexArray).toBe(0);
        expect(ctx.stats.resourceCounts.program).toBe(5);
        expect(ctx.stats.resourceCounts.shader).toBe(10);

        ctx.resources.destroy();
        expect(ctx.stats.resourceCounts.program).toBe(0);
        expect(ctx.stats.resourceCounts.shader).toBe(0);
    });

    it('transparency', async () => {
        const [width, height] = [32, 32];
        const gl = createGl(width, height, { preserveDrawingBuffer: true });
        const { ctx } = createRenderer(gl);
        const points = createPoints();

        const sceneBlended = Scene.create(ctx, 'blended');
        sceneBlended.add(points);
        sceneBlended.commit();

        const sceneWboit = Scene.create(ctx, 'wboit');
        sceneWboit.add(points);
        sceneWboit.commit();

        const sceneDpoit = Scene.create(ctx, 'dpoit');
        sceneDpoit.add(points);
        sceneDpoit.commit();

        expect(ctx.stats.resourceCounts.attribute).toBe(ctx.isWebGL2 ? 12 : 15);
        expect(ctx.stats.resourceCounts.texture).toBe(28);
        expect(ctx.stats.resourceCounts.vertexArray).toBe(ctx.extensions.vertexArrayObject ? 15 : 0);
        expect(ctx.stats.resourceCounts.program).toBe(7);
        expect(ctx.stats.resourceCounts.shader).toBe(14);
    });
});