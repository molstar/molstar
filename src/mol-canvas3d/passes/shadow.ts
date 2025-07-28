/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { TextureSpec, Values, UniformSpec, DefineSpec } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Texture } from '../../mol-gl/webgl/texture';
import { ValueCell } from '../../mol-util';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { createComputeRenderable, ComputeRenderable } from '../../mol-gl/renderable';
import { Mat4, Vec2, Vec3, Vec4 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { ICamera } from '../../mol-canvas3d/camera';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { isTimingMode } from '../../mol-util/debug';
import { Light } from '../../mol-gl/renderer';
import { shadows_frag } from '../../mol-gl/shader/shadows.frag';
import { PostprocessingProps } from './postprocessing';

export const ShadowParams = {
    steps: PD.Numeric(1, { min: 1, max: 64, step: 1 }),
    maxDistance: PD.Numeric(3, { min: 0, max: 256, step: 1 }),
    tolerance: PD.Numeric(1.0, { min: 0.0, max: 10.0, step: 0.1 }),
};

export type ShadowProps = PD.Values<typeof ShadowParams>

export class ShadowPass {
    static isEnabled(props: PostprocessingProps) {
        return props.enabled && props.shadow.name !== 'off';
    }

    readonly target: RenderTarget;
    private readonly renderable: ShadowsRenderable;

    constructor(readonly webgl: WebGLContext, width: number, height: number, depthTextureOpaque: Texture) {
        this.target = webgl.createRenderTarget(width, height, false);
        this.renderable = getShadowsRenderable(webgl, depthTextureOpaque);
    }

    setSize(width: number, height: number) {
        const [w, h] = this.renderable.values.uTexSize.ref.value;
        if (width !== w || height !== h) {
            this.target.setSize(width, height);
            ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
        }
    }

    update(camera: ICamera, light: Light, ambientColor: Vec3, props: ShadowProps) {
        let needsUpdateShadows = false;

        const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;

        const invProjection = Mat4.identity();
        Mat4.invert(invProjection, camera.projection);

        const [w, h] = this.renderable.values.uTexSize.ref.value;
        const v = camera.viewport;

        ValueCell.update(this.renderable.values.uProjection, camera.projection);
        ValueCell.update(this.renderable.values.uInvProjection, invProjection);

        Vec4.set(this.renderable.values.uBounds.ref.value,
            v.x / w,
            v.y / h,
            (v.x + v.width) / w,
            (v.y + v.height) / h
        );
        ValueCell.update(this.renderable.values.uBounds, this.renderable.values.uBounds.ref.value);

        ValueCell.updateIfChanged(this.renderable.values.uNear, camera.near);
        ValueCell.updateIfChanged(this.renderable.values.uFar, camera.far);
        if (this.renderable.values.dOrthographic.ref.value !== orthographic) {
            ValueCell.update(this.renderable.values.dOrthographic, orthographic);
            needsUpdateShadows = true;
        }

        ValueCell.updateIfChanged(this.renderable.values.uMaxDistance, props.maxDistance * camera.state.scale);
        ValueCell.updateIfChanged(this.renderable.values.uTolerance, props.tolerance * camera.state.scale);
        if (this.renderable.values.dSteps.ref.value !== props.steps) {
            ValueCell.update(this.renderable.values.dSteps, props.steps);
            needsUpdateShadows = true;
        }

        ValueCell.update(this.renderable.values.uLightDirection, light.direction);
        ValueCell.update(this.renderable.values.uLightColor, light.color);
        if (this.renderable.values.dLightCount.ref.value !== light.count) {
            ValueCell.update(this.renderable.values.dLightCount, light.count);
            needsUpdateShadows = true;
        }
        ValueCell.update(this.renderable.values.uAmbientColor, ambientColor);

        if (needsUpdateShadows) {
            this.renderable.update();
        }
    }

    render() {
        if (isTimingMode) this.webgl.timer.mark('ShadowPass.render');
        this.target.bind();
        this.renderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('ShadowPass.render');
    }
}

const ShadowsSchema = {
    ...QuadSchema,
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    uProjection: UniformSpec('m4'),
    uInvProjection: UniformSpec('m4'),
    uBounds: UniformSpec('v4'),

    dOrthographic: DefineSpec('number'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),

    dSteps: DefineSpec('number'),
    uMaxDistance: UniformSpec('f'),
    uTolerance: UniformSpec('f'),

    uLightDirection: UniformSpec('v3[]'),
    uLightColor: UniformSpec('v3[]'),
    dLightCount: DefineSpec('number'),
    uAmbientColor: UniformSpec('v3'),
};
type ShadowsRenderable = ComputeRenderable<Values<typeof ShadowsSchema>>

function getShadowsRenderable(ctx: WebGLContext, depthTexture: Texture): ShadowsRenderable {
    const width = depthTexture.getWidth();
    const height = depthTexture.getHeight();

    const values: Values<typeof ShadowsSchema> = {
        ...QuadValues,
        tDepth: ValueCell.create(depthTexture),
        uTexSize: ValueCell.create(Vec2.create(width, height)),

        uProjection: ValueCell.create(Mat4.identity()),
        uInvProjection: ValueCell.create(Mat4.identity()),
        uBounds: ValueCell.create(Vec4()),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),

        dSteps: ValueCell.create(1),
        uMaxDistance: ValueCell.create(3.0),
        uTolerance: ValueCell.create(1.0),

        uLightDirection: ValueCell.create([]),
        uLightColor: ValueCell.create([]),
        dLightCount: ValueCell.create(0),
        uAmbientColor: ValueCell.create(Vec3()),
    };

    const schema = { ...ShadowsSchema };
    const shaderCode = ShaderCode('shadows', quad_vert, shadows_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}