/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import Renderer from '../../mol-gl/renderer';
import Scene from '../../mol-gl/scene';
import { BoundingSphereHelper } from '../helper/bounding-sphere-helper';
import { Texture } from '../../mol-gl/webgl/texture';
import { ICamera } from '../camera';
import { CameraHelper, CameraHelperParams } from '../helper/camera-helper';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { HandleHelper } from '../helper/handle-helper';
import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { DefineSpec, TextureSpec, UniformSpec, Values } from '../../mol-gl/renderable/schema';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { ShaderCode } from '../../mol-gl/shader-code';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { ValueCell } from '../../mol-util';
import { Vec2 } from '../../mol-math/linear-algebra';
import { StereoCamera } from '../camera/stereo';

import quad_vert from '../../mol-gl/shader/quad.vert';
import depthMerge_frag from '../../mol-gl/shader/depth-merge.frag';

const DepthMergeSchema = {
    ...QuadSchema,
    tDepthPrimitives: TextureSpec('texture', 'depth', 'ushort', 'nearest'),
    tDepthVolumes: TextureSpec('texture', 'depth', 'ushort', 'nearest'),
    uTexSize: UniformSpec('v2'),
    dPackedDepth: DefineSpec('boolean'),
};
const DepthMergeShaderCode = ShaderCode('depth-merge', quad_vert, depthMerge_frag);
type DepthMergeRenderable = ComputeRenderable<Values<typeof DepthMergeSchema>>

function getDepthMergeRenderable(ctx: WebGLContext, depthTexturePrimitives: Texture, depthTextureVolumes: Texture, packedDepth: boolean): DepthMergeRenderable {
    const values: Values<typeof DepthMergeSchema> = {
        ...QuadValues,
        tDepthPrimitives: ValueCell.create(depthTexturePrimitives),
        tDepthVolumes: ValueCell.create(depthTextureVolumes),
        uTexSize: ValueCell.create(Vec2.create(depthTexturePrimitives.getWidth(), depthTexturePrimitives.getHeight())),
        dPackedDepth: ValueCell.create(packedDepth),
    };

    const schema = { ...DepthMergeSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', DepthMergeShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

export const DrawPassParams = {
    cameraHelper: PD.Group(CameraHelperParams)
};
export const DefaultDrawPassProps = PD.getDefaultValues(DrawPassParams);
export type DrawPassProps = PD.Values<typeof DrawPassParams>

export class DrawPass {
    readonly colorTarget: RenderTarget
    readonly depthTexture: Texture
    readonly depthTexturePrimitives: Texture
    readonly packedDepth: boolean

    readonly cameraHelper: CameraHelper

    private depthTarget: RenderTarget
    private depthTargetPrimitives: RenderTarget | null
    private depthTargetVolumes: RenderTarget | null
    private depthTextureVolumes: Texture
    private depthMerge: DepthMergeRenderable

    isStereo = false

    constructor(private webgl: WebGLContext, private renderer: Renderer, private scene: Scene, private camera: { standard: ICamera, stereo?: StereoCamera }, private debugHelper: BoundingSphereHelper, private handleHelper: HandleHelper, props: Partial<DrawPassProps> = {}) {
        const { extensions, resources } = webgl;
        const { width, height } = camera.standard.viewport;

        this.colorTarget = webgl.createRenderTarget(width, height);
        this.packedDepth = !extensions.depthTexture;

        this.depthTarget = webgl.createRenderTarget(width, height);
        this.depthTexture = this.depthTarget.texture;

        this.depthTargetPrimitives = this.packedDepth ? webgl.createRenderTarget(width, height) : null;
        this.depthTargetVolumes = this.packedDepth ? webgl.createRenderTarget(width, height) : null;

        this.depthTexturePrimitives = this.depthTargetPrimitives ? this.depthTargetPrimitives.texture : resources.texture('image-depth', 'depth', 'ushort', 'nearest');
        this.depthTextureVolumes = this.depthTargetVolumes ? this.depthTargetVolumes.texture : resources.texture('image-depth', 'depth', 'ushort', 'nearest');
        if (!this.packedDepth) {
            this.depthTexturePrimitives.define(width, height);
            this.depthTextureVolumes.define(width, height);
        }
        this.depthMerge = getDepthMergeRenderable(webgl, this.depthTexturePrimitives, this.depthTextureVolumes, this.packedDepth);

        const p = { ...DefaultDrawPassProps, ...props };
        this.cameraHelper = new CameraHelper(webgl, p.cameraHelper);
    }

    setSize(width: number, height: number) {
        this.colorTarget.setSize(width, height);
        this.depthTarget.setSize(width, height);

        if (this.depthTargetPrimitives) {
            this.depthTargetPrimitives.setSize(width, height);
        } else {
            this.depthTexturePrimitives.define(width, height);
        }

        if (this.depthTargetVolumes) {
            this.depthTargetVolumes.setSize(width, height);
        } else {
            this.depthTextureVolumes.define(width, height);
        }

        ValueCell.update(this.depthMerge.values.uTexSize, Vec2.set(this.depthMerge.values.uTexSize.ref.value, width, height));
    }

    setProps(props: Partial<DrawPassProps>) {
        if (props.cameraHelper) this.cameraHelper.setProps(props.cameraHelper);
    }

    get props(): DrawPassProps {
        return {
            cameraHelper: { ...this.cameraHelper.props }
        };
    }

    render(toDrawingBuffer: boolean, transparentBackground: boolean) {
        if (this.isStereo && this.camera.stereo) {
            this._render(this.camera.stereo.left, toDrawingBuffer, transparentBackground);
            this._render(this.camera.stereo.right, toDrawingBuffer, transparentBackground);
        } else {
            this._render(this.camera.standard, toDrawingBuffer, transparentBackground);
        }
    }

    private _render(camera: ICamera, toDrawingBuffer: boolean, transparentBackground: boolean) {
        const { x, y, width, height } = camera.viewport;
        if (toDrawingBuffer) {
            this.webgl.unbindFramebuffer();
            this.renderer.setViewport(x, y, width, height);
        } else {
            this.colorTarget.bind();
            this.renderer.setViewport(x, y, width, height);
            if (!this.packedDepth) {
                this.depthTexturePrimitives.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
            }
        }

        this.renderer.render(this.scene.primitives, camera, 'color', true, transparentBackground, null);

        // do a depth pass if not rendering to drawing buffer and
        // extensions.depthTexture is unsupported (i.e. depthTarget is set)
        if (!toDrawingBuffer && this.depthTargetPrimitives) {
            this.depthTargetPrimitives.bind();
            this.renderer.setViewport(x, y, width, height);
            this.renderer.render(this.scene.primitives, camera, 'depth', true, transparentBackground, null);
            this.colorTarget.bind();
        }

        // do direct-volume rendering
        if (!toDrawingBuffer) {
            if (!this.packedDepth) {
                this.depthTextureVolumes.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
                this.webgl.state.depthMask(true);
                this.webgl.gl.clear(this.webgl.gl.DEPTH_BUFFER_BIT);
            }
            this.renderer.render(this.scene.volumes, camera, 'color', false, transparentBackground, this.depthTexturePrimitives);

            // do volume depth pass if extensions.depthTexture is unsupported (i.e. depthTarget is set)
            if (this.depthTargetVolumes) {
                this.depthTargetVolumes.bind();
                this.renderer.render(this.scene.volumes, camera, 'depth', true, transparentBackground, this.depthTexturePrimitives);
                this.colorTarget.bind();
            }
        }

        // merge depths from primitive and volume rendering
        if (!toDrawingBuffer) {
            this.depthMerge.update();
            this.depthTarget.bind();
            this.renderer.setViewport(x, y, width, height);
            this.webgl.state.disable(this.webgl.gl.SCISSOR_TEST);
            this.webgl.state.disable(this.webgl.gl.BLEND);
            this.webgl.state.disable(this.webgl.gl.DEPTH_TEST);
            this.webgl.state.depthMask(false);
            this.webgl.state.clearColor(1, 1, 1, 1);
            this.webgl.gl.clear(this.webgl.gl.COLOR_BUFFER_BIT);
            this.depthMerge.render();
            this.colorTarget.bind();
            this.renderer.setViewport(x, y, width, height);
        }

        if (this.debugHelper.isEnabled) {
            this.debugHelper.syncVisibility();
            this.renderer.render(this.debugHelper.scene, camera, 'color', false, transparentBackground, null);
        }
        if (this.handleHelper.isEnabled) {
            this.renderer.render(this.handleHelper.scene, camera, 'color', false, transparentBackground, null);
        }
        if (this.cameraHelper.isEnabled) {
            this.cameraHelper.update(camera);
            this.renderer.render(this.cameraHelper.scene, this.cameraHelper.camera, 'color', false, transparentBackground, null);
        }
    }
}