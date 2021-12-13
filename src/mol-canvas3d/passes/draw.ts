/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { createNullRenderTarget, RenderTarget } from '../../mol-gl/webgl/render-target';
import { Renderer } from '../../mol-gl/renderer';
import { Scene } from '../../mol-gl/scene';
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

import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { depthMerge_frag } from '../../mol-gl/shader/depth-merge.frag';
import { StereoCamera } from '../camera/stereo';
import { WboitPass } from './wboit';
import { AntialiasingPass, PostprocessingPass, PostprocessingProps } from './postprocessing';
import { MarkingPass, MarkingProps } from './marking';
import { CopyRenderable, createCopyRenderable } from '../../mol-gl/compute/util';

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
    private readonly drawTarget: RenderTarget

    readonly colorTarget: RenderTarget
    readonly depthTexture: Texture
    readonly depthTexturePrimitives: Texture

    readonly packedDepth: boolean

    private depthTarget: RenderTarget
    private depthTargetPrimitives: RenderTarget | null
    private depthTargetVolumes: RenderTarget | null
    private depthTextureVolumes: Texture
    private depthMerge: DepthMergeRenderable

    private copyFboTarget: CopyRenderable
    private copyFboPostprocessing: CopyRenderable

    private wboit: WboitPass | undefined
    private readonly marking: MarkingPass
    readonly postprocessing: PostprocessingPass
    private readonly antialiasing: AntialiasingPass

    get wboitEnabled() {
        return !!this.wboit?.supported;
    }

    constructor(private webgl: WebGLContext, width: number, height: number, enableWboit: boolean) {
        const { extensions, resources, isWebGL2 } = webgl;

        this.drawTarget = createNullRenderTarget(webgl.gl);

        this.colorTarget = webgl.createRenderTarget(width, height, true, 'uint8', 'linear');
        this.packedDepth = !extensions.depthTexture;

        this.depthTarget = webgl.createRenderTarget(width, height);
        this.depthTexture = this.depthTarget.texture;

        this.depthTargetPrimitives = this.packedDepth ? webgl.createRenderTarget(width, height) : null;
        this.depthTargetVolumes = this.packedDepth ? webgl.createRenderTarget(width, height) : null;

        this.depthTexturePrimitives = this.depthTargetPrimitives ? this.depthTargetPrimitives.texture : resources.texture('image-depth', 'depth', isWebGL2 ? 'float' : 'ushort', 'nearest');
        this.depthTextureVolumes = this.depthTargetVolumes ? this.depthTargetVolumes.texture : resources.texture('image-depth', 'depth', isWebGL2 ? 'float' : 'ushort', 'nearest');
        if (!this.packedDepth) {
            this.depthTexturePrimitives.define(width, height);
            this.depthTextureVolumes.define(width, height);
        }
        this.depthMerge = getDepthMergeRenderable(webgl, this.depthTexturePrimitives, this.depthTextureVolumes, this.packedDepth);

        this.wboit = enableWboit ? new WboitPass(webgl, width, height) : undefined;
        this.marking = new MarkingPass(webgl, width, height);
        this.postprocessing = new PostprocessingPass(webgl, this);
        this.antialiasing = new AntialiasingPass(webgl, this);

        this.copyFboTarget = createCopyRenderable(webgl, this.colorTarget.texture);
        this.copyFboPostprocessing = createCopyRenderable(webgl, this.postprocessing.target.texture);
    }

    reset() {
        this.wboit?.reset();
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

            if (this.depthTargetVolumes) {
                this.depthTargetVolumes.setSize(width, height);
            } else {
                this.depthTextureVolumes.define(width, height);
            }

            ValueCell.update(this.depthMerge.values.uTexSize, Vec2.set(this.depthMerge.values.uTexSize.ref.value, width, height));

            ValueCell.update(this.copyFboTarget.values.uTexSize, Vec2.set(this.copyFboTarget.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.copyFboPostprocessing.values.uTexSize, Vec2.set(this.copyFboPostprocessing.values.uTexSize.ref.value, width, height));

            if (this.wboit?.supported) {
                this.wboit.setSize(width, height);
            }

            this.marking.setSize(width, height);
            this.postprocessing.setSize(width, height);
            this.antialiasing.setSize(width, height);
        }
    }

    private _depthMerge() {
        const { state, gl } = this.webgl;

        this.depthMerge.update();
        this.depthTarget.bind();
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.disable(gl.CULL_FACE);
        state.depthMask(false);
        state.clearColor(1, 1, 1, 1);
        gl.clear(gl.COLOR_BUFFER_BIT);
        this.depthMerge.render();
    }

    private _renderWboit(renderer: Renderer, camera: ICamera, scene: Scene, transparentBackground: boolean, postprocessingProps: PostprocessingProps) {
        if (!this.wboit?.supported) throw new Error('expected wboit to be supported');

        this.colorTarget.bind();
        renderer.clear(true);

        // render opaque primitives
        this.depthTexturePrimitives.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
        this.colorTarget.bind();
        renderer.clearDepth();
        renderer.renderWboitOpaque(scene.primitives, camera, null);

        // render opaque volumes
        this.depthTextureVolumes.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
        this.colorTarget.bind();
        renderer.clearDepth();
        renderer.renderWboitOpaque(scene.volumes, camera, this.depthTexturePrimitives);

        // merge depth of opaque primitives and volumes
        this._depthMerge();

        if (PostprocessingPass.isEnabled(postprocessingProps)) {
            this.postprocessing.render(camera, false, transparentBackground, renderer.props.backgroundColor, postprocessingProps);
        }

        // render transparent primitives and volumes
        this.wboit.bind();
        renderer.renderWboitTransparent(scene.primitives, camera, this.depthTexture);
        renderer.renderWboitTransparent(scene.volumes, camera, this.depthTexture);

        // evaluate wboit
        if (PostprocessingPass.isEnabled(postprocessingProps)) {
            this.depthTexturePrimitives.attachFramebuffer(this.postprocessing.target.framebuffer, 'depth');
            this.postprocessing.target.bind();
        } else {
            this.depthTexturePrimitives.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
            this.colorTarget.bind();
        }
        this.wboit.render();
    }

    private _renderBlended(renderer: Renderer, camera: ICamera, scene: Scene, toDrawingBuffer: boolean, transparentBackground: boolean, postprocessingProps: PostprocessingProps) {
        if (toDrawingBuffer) {
            this.drawTarget.bind();
        } else {
            this.colorTarget.bind();
            if (!this.packedDepth) {
                this.depthTexturePrimitives.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
            }
        }

        renderer.clear(true);
        renderer.renderBlendedOpaque(scene.primitives, camera, null);

        if (!toDrawingBuffer) {
            // do a depth pass if not rendering to drawing buffer and
            // extensions.depthTexture is unsupported (i.e. depthTarget is set)
            if (this.depthTargetPrimitives) {
                this.depthTargetPrimitives.bind();
                renderer.clear(false);
                // TODO: this should only render opaque
                renderer.renderDepth(scene.primitives, camera, null);
                this.colorTarget.bind();
            }

            // do direct-volume rendering
            if (!this.packedDepth) {
                this.depthTextureVolumes.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
                renderer.clearDepth(); // from previous frame
            }
            renderer.renderBlendedVolumeOpaque(scene.volumes, camera, this.depthTexturePrimitives);

            // do volume depth pass if extensions.depthTexture is unsupported (i.e. depthTarget is set)
            if (this.depthTargetVolumes) {
                this.depthTargetVolumes.bind();
                renderer.clear(false);
                renderer.renderDepth(scene.volumes, camera, this.depthTexturePrimitives);
                this.colorTarget.bind();
            }

            // merge depths from primitive and volume rendering
            this._depthMerge();
            this.colorTarget.bind();

            if (PostprocessingPass.isEnabled(postprocessingProps)) {
                this.postprocessing.render(camera, false, transparentBackground, renderer.props.backgroundColor, postprocessingProps);
            }
            renderer.renderBlendedVolumeTransparent(scene.volumes, camera, this.depthTexturePrimitives);

            const target = PostprocessingPass.isEnabled(postprocessingProps)
                ? this.postprocessing.target : this.colorTarget;
            if (!this.packedDepth) {
                this.depthTexturePrimitives.attachFramebuffer(target.framebuffer, 'depth');
            }
            target.bind();
        }

        renderer.renderBlendedTransparent(scene.primitives, camera, null);
    }

    private _render(renderer: Renderer, camera: ICamera, scene: Scene, helper: Helper, toDrawingBuffer: boolean, transparentBackground: boolean, postprocessingProps: PostprocessingProps, markingProps: MarkingProps) {
        const volumeRendering = scene.volumes.renderables.length > 0;
        const postprocessingEnabled = PostprocessingPass.isEnabled(postprocessingProps);
        const antialiasingEnabled = AntialiasingPass.isEnabled(postprocessingProps);
        const markingEnabled = MarkingPass.isEnabled(markingProps);

        const { x, y, width, height } = camera.viewport;
        renderer.setViewport(x, y, width, height);
        renderer.update(camera);

        if (transparentBackground && !antialiasingEnabled && toDrawingBuffer) {
            this.drawTarget.bind();
            renderer.clear(false);
        }

        if (this.wboitEnabled) {
            this._renderWboit(renderer, camera, scene, transparentBackground, postprocessingProps);
        } else {
            this._renderBlended(renderer, camera, scene, !volumeRendering && !postprocessingEnabled && !antialiasingEnabled && toDrawingBuffer, transparentBackground, postprocessingProps);
        }

        if (postprocessingEnabled) {
            this.postprocessing.target.bind();
        } else if (!toDrawingBuffer || volumeRendering || this.wboitEnabled) {
            this.colorTarget.bind();
        } else {
            this.drawTarget.bind();
        }

        if (markingEnabled) {
            const markingDepthTest = markingProps.ghostEdgeStrength < 1;
            if (markingDepthTest) {
                this.marking.depthTarget.bind();
                renderer.clear(false);
                renderer.renderMarkingDepth(scene.primitives, camera, null);
            }

            this.marking.maskTarget.bind();
            renderer.clear(false);
            renderer.renderMarkingMask(scene.primitives, camera, markingDepthTest ? this.marking.depthTarget.texture : null);

            this.marking.update(markingProps);
            this.marking.render(camera.viewport, postprocessingEnabled ? this.postprocessing.target : this.colorTarget);
        }

        if (helper.debug.isEnabled) {
            helper.debug.syncVisibility();
            renderer.renderBlended(helper.debug.scene, camera, null);
        }
        if (helper.handle.isEnabled) {
            renderer.renderBlended(helper.handle.scene, camera, null);
        }
        if (helper.camera.isEnabled) {
            helper.camera.update(camera);
            renderer.update(helper.camera.camera);
            renderer.renderBlended(helper.camera.scene, helper.camera.camera, null);
        }

        if (antialiasingEnabled) {
            this.antialiasing.render(camera, toDrawingBuffer, postprocessingProps);
        } else if (toDrawingBuffer) {
            this.drawTarget.bind();

            this.webgl.state.disable(this.webgl.gl.DEPTH_TEST);
            if (postprocessingEnabled) {
                this.copyFboPostprocessing.render();
            } else if (volumeRendering || this.wboitEnabled) {
                this.copyFboTarget.render();
            }
        }

        this.webgl.gl.flush();
    }

    render(renderer: Renderer, camera: Camera | StereoCamera, scene: Scene, helper: Helper, toDrawingBuffer: boolean, transparentBackground: boolean, postprocessingProps: PostprocessingProps, markingProps: MarkingProps) {
        renderer.setTransparentBackground(transparentBackground);
        renderer.setDrawingBufferSize(this.colorTarget.getWidth(), this.colorTarget.getHeight());
        renderer.setPixelRatio(this.webgl.pixelRatio);

        if (StereoCamera.is(camera)) {
            this._render(renderer, camera.left, scene, helper, toDrawingBuffer, transparentBackground, postprocessingProps, markingProps);
            this._render(renderer, camera.right, scene, helper, toDrawingBuffer, transparentBackground, postprocessingProps, markingProps);
        } else {
            this._render(renderer, camera, scene, helper, toDrawingBuffer, transparentBackground, postprocessingProps, markingProps);
        }
    }

    getColorTarget(postprocessingProps: PostprocessingProps): RenderTarget {
        if (AntialiasingPass.isEnabled(postprocessingProps)) {
            return this.antialiasing.target;
        } else if (PostprocessingPass.isEnabled(postprocessingProps)) {
            return this.postprocessing.target;
        }
        return this.colorTarget;
    }
}