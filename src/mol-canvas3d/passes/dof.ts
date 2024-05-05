/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { TextureSpec, UniformSpec, DefineSpec, Values } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { Texture } from '../../mol-gl/webgl/texture';
import { Mat4, Vec2, Vec3, Vec4 } from '../../mol-math/linear-algebra';
import { ValueCell } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { dof_frag } from '../../mol-gl/shader/dof.frag';
import { Viewport } from '../camera/util';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { isTimingMode } from '../../mol-util/debug';
import { ICamera } from '../../mol-canvas3d/camera';
import { Sphere3D } from '../../mol-math/geometry';

export const DofParams = {
    blurSize: PD.Numeric(9, { min: 1, max: 32, step: 1 }),
    blurSpread: PD.Numeric(1.0, { min: 0.0, max: 10.0, step: 0.1 }),
    inFocus: PD.Numeric(0.0, { min: -5000.0, max: 5000.0, step: 1.0 }, { description: 'Distance from the scene center that will be in focus' }),
    PPM: PD.Numeric(20.0, { min: 0.0, max: 5000.0, step: 0.1 }, { description: 'Size of the area that will be in focus' }),
    center: PD.Select('scene-center', PD.arrayToOptions(['scene-center', 'camera-target'])),
    mode: PD.Select('plane', PD.arrayToOptions(['plane', 'sphere'])),
};

export type DofProps = PD.Values<typeof DofParams>

export class DofPass {
    private readonly renderable: DofRenderable;

    constructor(private webgl: WebGLContext, input: Texture, depth: Texture) {
        this.renderable = getDofRenderable(webgl, input, depth);
    }

    private updateState(viewport: Viewport) {
        const { gl, state } = this.webgl;

        state.enable(gl.SCISSOR_TEST);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);

        const { x, y, width, height } = viewport;
        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);

        state.clearColor(0, 0, 0, 1);
        gl.clear(gl.COLOR_BUFFER_BIT);
    }

    setSize(width: number, height: number) {
        ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
    }

    update(camera: ICamera, input: Texture, depth: Texture, props: DofProps, sphere: Sphere3D) {
        let needsUpdate = false;
        if (this.renderable.values.tColor.ref.value !== input) {
            ValueCell.update(this.renderable.values.tColor, input);
            needsUpdate = true;
        }
        if (this.renderable.values.tDepth.ref.value !== depth) {
            ValueCell.update(this.renderable.values.tDepth, depth);
            needsUpdate = true;
        }
        const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;
        const invProjection = Mat4.identity();
        Mat4.invert(invProjection, camera.projection);
        const wolrd_center = (props.center === 'scene-center' ? sphere.center : camera.state.target);
        const [w, h] = this.renderable.values.uTexSize.ref.value;
        const v = camera.viewport;
        const distance = Vec3.distance(camera.state.position, wolrd_center);
        const inFocus = distance + props.inFocus;
        // transform  center in view space
        let center = Vec3();
        center = Vec3.transformMat4(center, wolrd_center, camera.view);
        ValueCell.update(this.renderable.values.uProjection, camera.projection);
        ValueCell.update(this.renderable.values.uInvProjection, invProjection);
        ValueCell.update(this.renderable.values.uMode, props.mode === 'sphere' ? 1 : 0);

        Vec4.set(this.renderable.values.uBounds.ref.value,
            v.x / w,
            v.y / h,
            (v.x + v.width) / w,
            (v.y + v.height) / h
        );
        ValueCell.update(this.renderable.values.uBounds, this.renderable.values.uBounds.ref.value);

        ValueCell.updateIfChanged(this.renderable.values.uNear, camera.near);
        ValueCell.updateIfChanged(this.renderable.values.uFar, camera.far);
        ValueCell.updateIfChanged(this.renderable.values.dOrthographic, orthographic);

        if (this.renderable.values.uCenter.ref.value !== center) needsUpdate = true;
        ValueCell.update(this.renderable.values.uCenter, center);

        if (this.renderable.values.blurSize.ref.value !== props.blurSize) needsUpdate = true;
        ValueCell.update(this.renderable.values.blurSize, props.blurSize);
        if (this.renderable.values.blurSpread.ref.value !== props.blurSpread) needsUpdate = true;
        ValueCell.update(this.renderable.values.blurSpread, props.blurSpread);
        if (this.renderable.values.inFocus.ref.value !== inFocus) needsUpdate = true;
        ValueCell.update(this.renderable.values.inFocus, inFocus);
        if (this.renderable.values.PPM.ref.value !== props.PPM) needsUpdate = true;
        ValueCell.update(this.renderable.values.PPM, props.PPM);

        this.renderable.update();
        ValueCell.update(this.renderable.values.PPM, props.PPM);
        ValueCell.update(this.renderable.values.inFocus, inFocus);
        ValueCell.update(this.renderable.values.blurSpread, props.blurSpread);
        ValueCell.update(this.renderable.values.blurSize, props.blurSize);

        if (needsUpdate) {
            this.renderable.update();
        }
    }

    render(viewport: Viewport, target: RenderTarget | undefined) {
        if (isTimingMode) this.webgl.timer.mark('DofPass.render');
        if (target) {
            target.bind();
        } else {
            this.webgl.unbindFramebuffer();
        }
        this.updateState(viewport);
        this.renderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('DofPass.render');
    }
}

//

const DofSchema = {
    ...QuadSchema,
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    uProjection: UniformSpec('m4'),
    uInvProjection: UniformSpec('m4'),
    uBounds: UniformSpec('v4'),
    uCenter: UniformSpec('v3'),
    uMode: UniformSpec('i'),

    dOrthographic: DefineSpec('number'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),

    inFocus: UniformSpec('f'),
    PPM: UniformSpec('f'),
    blurSize: UniformSpec('i'),
    blurSpread: UniformSpec('f'),
};

const DofShaderCode = ShaderCode('dof', quad_vert, dof_frag);
type DofRenderable = ComputeRenderable<Values<typeof DofSchema>>

function getDofRenderable(ctx: WebGLContext, colorTexture: Texture, depthTexture: Texture): DofRenderable {
    const width = colorTexture.getWidth();
    const height = colorTexture.getHeight();

    const values: Values<typeof DofSchema> = {
        ...QuadValues,
        tDepth: ValueCell.create(depthTexture),
        tColor: ValueCell.create(colorTexture),
        uTexSize: ValueCell.create(Vec2.create(width, height)),

        uProjection: ValueCell.create(Mat4.identity()),
        uInvProjection: ValueCell.create(Mat4.identity()),
        uBounds: ValueCell.create(Vec4()),
        uCenter: ValueCell.create(Vec3()),
        uMode: ValueCell.create(0),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),

        blurSize: ValueCell.create(5),
        blurSpread: ValueCell.create(300.0),
        inFocus: ValueCell.create(20.0),
        PPM: ValueCell.create(20.0),
    };

    const schema = { ...DofSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', DofShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}