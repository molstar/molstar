/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createGl } from './gl.shim';

import { Camera } from '../../mol-canvas3d/camera';
import { Vec3, Mat4, Vec4 } from '../../mol-math/linear-algebra';
import { ValueCell } from '../../mol-util';

import Renderer from '../renderer';
import { createValueColor } from '../../mol-geo/geometry/color-data';
import { createValueSize } from '../../mol-geo/geometry/size-data';
import { createContext } from '../webgl/context';
import { RenderableState } from '../renderable';
import { createRenderObject } from '../render-object';
import { PointsValues } from '../renderable/points';
import Scene from '../scene';
import { createEmptyMarkers } from '../../mol-geo/geometry/marker-data';
import { fillSerial } from '../../mol-util/array';
import { Color } from '../../mol-util/color';
import { Sphere3D } from '../../mol-math/geometry';
import { createEmptyOverpaint } from '../../mol-geo/geometry/overpaint-data';
import { createEmptyTransparency } from '../../mol-geo/geometry/transparency-data';
import { createEmptyClipping } from '../../mol-geo/geometry/clipping-data';

function createRenderer(gl: WebGLRenderingContext) {
    const ctx = createContext(gl);
    const camera = new Camera({
        position: Vec3.create(0, 0, 50)
    });
    const renderer = Renderer.create(ctx);
    return { ctx, camera, renderer };
}

function createPoints() {
    const aPosition = ValueCell.create(new Float32Array([0, -1, 0, -1, 0, 0, 1, 1, 0]));
    const aGroup = ValueCell.create(fillSerial(new Float32Array(3)));
    const aInstance = ValueCell.create(fillSerial(new Float32Array(1)));
    const color = createValueColor(Color(0xFF0000));
    const size = createValueSize(1);
    const marker = createEmptyMarkers();
    const overpaint = createEmptyOverpaint();
    const transparency = createEmptyTransparency();
    const clipping = createEmptyClipping();

    const aTransform = ValueCell.create(new Float32Array(16));
    const m4 = Mat4.identity();
    Mat4.toArray(m4, aTransform.ref.value, 0);
    const transform = ValueCell.create(new Float32Array(aTransform.ref.value));
    const extraTransform = ValueCell.create(new Float32Array(aTransform.ref.value));

    const boundingSphere = ValueCell.create(Sphere3D.create(Vec3.zero(), 2));
    const invariantBoundingSphere = ValueCell.create(Sphere3D.create(Vec3.zero(), 2));

    const values: PointsValues = {
        aPosition,
        aGroup,
        aTransform,
        aInstance,
        ...color,
        ...marker,
        ...size,
        ...overpaint,
        ...transparency,
        ...clipping,

        uAlpha: ValueCell.create(1.0),
        uInstanceCount: ValueCell.create(1),
        uGroupCount: ValueCell.create(3),
        uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere.ref.value)),

        alpha: ValueCell.create(1.0),
        drawCount: ValueCell.create(3),
        instanceCount: ValueCell.create(1),
        matrix: ValueCell.create(m4),
        transform,
        extraTransform,
        boundingSphere,
        invariantBoundingSphere,

        uSizeFactor: ValueCell.create(1),
        dPointSizeAttenuation: ValueCell.create(true),
        dPointFilledCircle: ValueCell.create(false),
        uPointEdgeBleach: ValueCell.create(0.5),
    };
    const state: RenderableState = {
        visible: true,
        alphaFactor: 1,
        pickable: true,
        opaque: true,
        writeDepth: true
    };

    return createRenderObject('points', values, state, -1);
}

describe('renderer', () => {
    it('basic', () => {
        const [ width, height ] = [ 32, 32 ];
        const gl = createGl(width, height, { preserveDrawingBuffer: true });
        const { ctx, renderer } = createRenderer(gl);

        expect(ctx.gl.canvas.width).toBe(32);
        expect(ctx.gl.canvas.height).toBe(32);

        expect(ctx.stats.resourceCounts.attribute).toBe(0);
        expect(ctx.stats.resourceCounts.texture).toBe(0);
        expect(ctx.stats.resourceCounts.vertexArray).toBe(0);
        expect(ctx.stats.resourceCounts.program).toBe(0);
        expect(ctx.stats.resourceCounts.shader).toBe(0);

        renderer.setViewport(0, 0, 64, 48);
        expect(ctx.gl.getParameter(ctx.gl.VIEWPORT)[2]).toBe(64);
        expect(ctx.gl.getParameter(ctx.gl.VIEWPORT)[3]).toBe(48);
    });

    it('points', async () => {
        const [ width, height ] = [ 32, 32 ];
        const gl = createGl(width, height, { preserveDrawingBuffer: true });
        const { ctx } = createRenderer(gl);
        const scene = Scene.create(ctx);

        const points = createPoints();

        scene.add(points);
        scene.commit();
        expect(ctx.stats.resourceCounts.attribute).toBe(4);
        expect(ctx.stats.resourceCounts.texture).toBe(6);
        expect(ctx.stats.resourceCounts.vertexArray).toBe(5);
        expect(ctx.stats.resourceCounts.program).toBe(5);
        expect(ctx.stats.resourceCounts.shader).toBe(10);

        scene.remove(points);
        scene.commit();
        expect(ctx.stats.resourceCounts.attribute).toBe(0);
        expect(ctx.stats.resourceCounts.texture).toBe(0);
        expect(ctx.stats.resourceCounts.vertexArray).toBe(0);
        expect(ctx.stats.resourceCounts.program).toBe(5);
        expect(ctx.stats.resourceCounts.shader).toBe(10);

        ctx.resources.destroy();
        expect(ctx.stats.resourceCounts.program).toBe(0);
        expect(ctx.stats.resourceCounts.shader).toBe(0);
    });
});