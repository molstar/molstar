/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import Renderer from '../../mol-gl/renderer';
import Scene from '../../mol-gl/scene';
import { Texture } from '../../mol-gl/webgl/texture';
import { Camera, ICamera } from '../camera';
import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { DefineSpec, TextureSpec, UniformSpec, Values } from '../../mol-gl/renderable/schema';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { ShaderCode } from '../../mol-gl/shader-code';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { ValueCell } from '../../mol-util';
import { Vec2 } from '../../mol-math/linear-algebra';
import { Helper } from '../helper/helper';

import quad_vert from '../../mol-gl/shader/quad.vert';
import depthMerge_frag from '../../mol-gl/shader/depth-merge.frag';
import { StereoCamera } from '../camera/stereo';

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

export class DrawPass {
    readonly colorTarget: RenderTarget
    readonly depthTexture: Texture
    readonly depthTexturePrimitives: Texture
    readonly depthTexturePrimitivesTransparent: Texture

    private readonly packedDepth: boolean
    private depthTarget: RenderTarget
    private depthTargetPrimitives: RenderTarget | null
    private depthTargetPrimitivesTransparent: RenderTarget | null
    private depthTargetVolumes: RenderTarget | null
    private depthTextureVolumes: Texture
    private depthMerge: DepthMergeRenderable

    constructor(private webgl: WebGLContext, width: number, height: number) {
        const { extensions, resources } = webgl;

        this.colorTarget = webgl.createRenderTarget(width, height);
        this.packedDepth = !extensions.depthTexture;

        this.depthTarget = webgl.createRenderTarget(width, height);
        this.depthTexture = this.depthTarget.texture;

        this.depthTargetPrimitives = this.packedDepth ? webgl.createRenderTarget(width, height) : null;
        this.depthTargetPrimitivesTransparent = this.packedDepth ? webgl.createRenderTarget(width, height) : null;
        this.depthTargetVolumes = this.packedDepth ? webgl.createRenderTarget(width, height) : null;

        this.depthTexturePrimitives = this.depthTargetPrimitives ? this.depthTargetPrimitives.texture : resources.texture('image-depth', 'depth', 'ushort', 'nearest');
        this.depthTexturePrimitivesTransparent = this.depthTargetPrimitivesTransparent ? this.depthTargetPrimitivesTransparent.texture : resources.texture('image-depth', 'depth', 'ushort', 'nearest');
        this.depthTextureVolumes = this.depthTargetVolumes ? this.depthTargetVolumes.texture : resources.texture('image-depth', 'depth', 'ushort', 'nearest');
        if (!this.packedDepth) {
            this.depthTexturePrimitives.define(width, height);
            this.depthTextureVolumes.define(width, height);
        }
        this.depthMerge = getDepthMergeRenderable(webgl, this.depthTexturePrimitives, this.depthTextureVolumes, this.packedDepth);
    }

    setSize(width: number, height: number) {
        const w = this.colorTarget.getWidth();
        const h = this.colorTarget.getHeight();

        if (width !== w || height !== h) {
            this.colorTarget.setSize(width, height);
            this.depthTarget.setSize(width, height);

            if (this.depthTargetPrimitives) {
                this.depthTargetPrimitives.setSize(width, height);
            } else {
                this.depthTexturePrimitives.define(width, height);
            }
            if (this.depthTargetPrimitivesTransparent) {
                this.depthTargetPrimitivesTransparent.setSize(width, height);
            } else {
                this.depthTexturePrimitivesTransparent.define(width, height);
            }

            if (this.depthTargetVolumes) {
                this.depthTargetVolumes.setSize(width, height);
            } else {
                this.depthTextureVolumes.define(width, height);
            }

            ValueCell.update(this.depthMerge.values.uTexSize, Vec2.set(this.depthMerge.values.uTexSize.ref.value, width, height));
        }
    }

    _render(renderer: Renderer, camera: ICamera, scene: Scene, helper: Helper, toDrawingBuffer: boolean, transparentBackground: boolean) {
        const { x, y, width, height } = camera.viewport;
        renderer.setViewport(x, y, width, height);

        let renderTarget;
        if (toDrawingBuffer) {
            renderTarget = null;
        } else {
            if (!this.packedDepth) {
                this.depthTexturePrimitives.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
            }
            renderTarget = this.colorTarget;
        }

        // do a depth pass if not rendering to drawing buffer and
        // extensions.depthTexture is unsupported (i.e. depthTarget is set)
        if (!toDrawingBuffer && this.depthTargetPrimitives) {
            renderer.render(this.depthTargetPrimitives, scene.primitives, camera, 'depth', true, transparentBackground, 1, null, false);
        }

        this.depthTexturePrimitives.attachFramebuffer(renderTarget!.framebuffer, 'depth');
        renderer.setViewport(0, 0, this.colorTarget.getWidth(), this.colorTarget.getHeight());
        renderer.render(renderTarget, scene.primitives, camera, 'color', true, transparentBackground, 1, null, false);

        if (helper.debug.isEnabled) {
            helper.debug.syncVisibility();
            renderer.render(renderTarget, helper.debug.scene, camera, 'color', false, transparentBackground, 1, null, false);
        }
        if (helper.handle.isEnabled) {
            renderer.render(renderTarget, helper.handle.scene, camera, 'color', false, transparentBackground, 1, null, false);
        }
        if (helper.camera.isEnabled) {
            helper.camera.update(camera);
            renderer.render(renderTarget, helper.camera.scene, helper.camera.camera, 'color', false, transparentBackground, 1, null, false);
        }

        if (!toDrawingBuffer && this.depthTargetPrimitivesTransparent) {
            renderer.render(this.depthTargetPrimitivesTransparent, scene.primitives, camera, 'depth', true, transparentBackground, 1, null, false);
        }

        if (renderTarget !== null) {
            // this.depthTexturePrimitivesTransparent.attachFramebuffer(renderTarget.framebuffer, 'depth');
            renderer.setViewport(0, 0, this.colorTarget.getWidth(), this.colorTarget.getHeight());
            renderer.render(renderTarget, scene.primitives, camera, 'color', false, transparentBackground, 1, this.depthTexturePrimitives, true);
        }

        if (helper.debug.isEnabled) {
            helper.debug.syncVisibility();
            renderer.render(renderTarget, helper.debug.scene, camera, 'color', false, transparentBackground, 1, this.depthTexturePrimitives, true);
        }
        if (helper.handle.isEnabled) {
            renderer.render(renderTarget, helper.handle.scene, camera, 'color', false, transparentBackground, 1, this.depthTexturePrimitives, true);
        }
        if (helper.camera.isEnabled) {
            helper.camera.update(camera);
            renderer.render(renderTarget, helper.camera.scene, helper.camera.camera, 'color', false, transparentBackground, 1, this.depthTexturePrimitives, true);
        }

        // do direct-volume rendering
        if (!toDrawingBuffer) {
            if (!this.packedDepth) {
                this.depthTextureVolumes.attachFramebuffer(renderTarget!.framebuffer, 'depth');
                this.webgl.state.depthMask(true);
                this.webgl.gl.viewport(x, y, width, height);
                this.webgl.gl.scissor(x, y, width, height);
                this.webgl.gl.clear(this.webgl.gl.DEPTH_BUFFER_BIT);
            }
            renderer.render(renderTarget, scene.volumes, camera, 'color', false, transparentBackground, 1, this.depthTexturePrimitives, true);

            // do volume depth pass if extensions.depthTexture is unsupported (i.e. depthTarget is set)
            if (this.depthTargetVolumes) {
                renderer.render(this.depthTargetVolumes, scene.volumes, camera, 'depth', true, transparentBackground, 1, this.depthTexturePrimitives, false);
            }
        }

        // merge depths from primitive and volume rendering
        if (!toDrawingBuffer) {
            this.depthMerge.update();
            this.depthTarget.bind();
            // this.webgl.state.disable(this.webgl.gl.SCISSOR_TEST);
            this.webgl.state.disable(this.webgl.gl.BLEND);
            this.webgl.state.disable(this.webgl.gl.DEPTH_TEST);
            this.webgl.state.disable(this.webgl.gl.CULL_FACE);
            this.webgl.state.depthMask(false);
            this.webgl.state.clearColor(1, 1, 1, 1);
            this.webgl.gl.viewport(x, y, width, height);
            this.webgl.gl.scissor(x, y, width, height);
            this.webgl.gl.clear(this.webgl.gl.COLOR_BUFFER_BIT);
            this.depthMerge.render();
            this.colorTarget.bind();
        }
    }

    render(renderer: Renderer, camera: Camera | StereoCamera, scene: Scene, helper: Helper, toDrawingBuffer: boolean, transparentBackground: boolean) {
        if (StereoCamera.is(camera)) {
            this._render(renderer, camera.left, scene, helper, toDrawingBuffer, transparentBackground);
            this._render(renderer, camera.right, scene, helper, toDrawingBuffer, transparentBackground);
        } else {
            this._render(renderer, camera, scene, helper, toDrawingBuffer, transparentBackground);
        }
    }
}