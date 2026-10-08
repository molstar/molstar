/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { Renderer } from '../../mol-gl/renderer';
import { Scene } from '../../mol-gl/scene';
import { Texture } from '../../mol-gl/webgl/texture';
import { ICamera } from '../camera';
import { Viewport } from '../camera/util';
import { ValueCell } from '../../mol-util';
import { Vec2 } from '../../mol-math/linear-algebra';

import { StereoCamera } from '../camera/stereo';
import { WboitPass } from './wboit';
import { DpoitPass } from './dpoit';
import { AntialiasingPass, PostprocessingPass, PostprocessingProps } from './postprocessing';
import { MarkingPass, MarkingProps, MarkingShading, SingleSample } from './marking';
import { CopyRenderable, createCopyRenderable } from '../../mol-gl/compute/util';
import { isDebugMode, isTimingMode } from '../../mol-util/debug';
import { AssetManager } from '../../mol-util/assets';
import { DofPass } from './dof';
import { BloomPass } from './bloom';
import { SsaoPass } from './ssao';
import { RenderContext } from '../util';

type Props = {
    postprocessing: PostprocessingProps;
    marking: MarkingProps;
    transparentBackground: boolean;
    dpoitIterations: number;
}

type TransparencyMode = 'wboit' | 'dpoit' | 'blended'

export class DrawPass {
    private readonly drawTarget: RenderTarget;

    readonly colorTarget: RenderTarget;
    readonly transparentColorTarget: RenderTarget;
    readonly depthTextureTransparent: Texture;
    readonly depthTextureOpaque: Texture;

    readonly packedDepth: boolean;

    readonly depthTargetTransparent: RenderTarget;
    private depthTargetOpaque: RenderTarget | null;

    private copyFbo: CopyRenderable;

    readonly wboit: WboitPass;
    readonly dpoit: DpoitPass;
    readonly marking: MarkingPass;
    readonly postprocessing: PostprocessingPass;
    private readonly ssaoShading: MarkingShading;
    readonly antialiasing: AntialiasingPass;
    readonly dof: DofPass;

    private transparencyMode: TransparencyMode = 'blended';
    setTransparency(transparency: 'wboit' | 'dpoit' | 'blended') {
        if (transparency === 'wboit') {
            this.transparencyMode = this.wboit.supported ? 'wboit' : 'blended';
            if (isDebugMode && !this.wboit.supported) {
                console.log('Missing "wboit" support, falling back to "blended".');
            }
        } else if (transparency === 'dpoit') {
            this.transparencyMode = this.dpoit.supported ? 'dpoit' : 'blended';
            if (isDebugMode && !this.dpoit.supported) {
                console.log('Missing "dpoit" support, falling back to "blended".');
            }
        } else {
            this.transparencyMode = 'blended';
        }
        this.depthTextureOpaque.detachFramebuffer(this.postprocessing.target.framebuffer, 'depth');
        this.marking.invalidate();
    }
    get transparency() {
        return this.transparencyMode;
    }

    constructor(private webgl: WebGLContext, assetManager: AssetManager, width: number, height: number, transparency: 'wboit' | 'dpoit' | 'blended') {
        const { extensions, resources, isWebGL2 } = webgl;
        this.drawTarget = webgl.createDrawTarget();
        this.colorTarget = webgl.createRenderTarget(width, height, 'depth-stencil', 'uint8', 'linear');
        this.transparentColorTarget = webgl.createRenderTarget(width, height, 'none', 'uint8', 'linear');

        this.packedDepth = !extensions.depthTexture;

        this.depthTargetTransparent = webgl.createRenderTarget(width, height, 'depth-stencil', 'uint8', 'nearest');
        this.depthTextureTransparent = this.depthTargetTransparent.texture;

        this.depthTargetOpaque = this.packedDepth ? webgl.createRenderTarget(width, height, 'depth-stencil') : null;

        this.depthTextureOpaque = this.depthTargetOpaque ? this.depthTargetOpaque.texture : resources.texture('image-depth', 'depth-stencil', isWebGL2 ? 'float-stencil' : 'uint24-8', 'nearest');
        if (!this.packedDepth) {
            this.depthTextureOpaque.define(width, height);
        }

        this.wboit = new WboitPass(webgl, width, height);
        this.dpoit = new DpoitPass(webgl, width, height, this.colorTarget.depthRenderbuffer);
        this.marking = new MarkingPass(webgl, width, height);
        this.postprocessing = new PostprocessingPass(webgl, assetManager, this);
        this.ssaoShading = { name: 'ssao', ssao: this.postprocessing.ssao.ssaoDepthTexture };
        this.antialiasing = new AntialiasingPass(webgl, width, height);
        this.dof = new DofPass(webgl, width, height);

        this.copyFbo = createCopyRenderable(webgl, this.colorTarget.texture);

        this.setTransparency(transparency);
    }

    getByteCount() {
        return (
            this.drawTarget.getByteCount() +
            this.colorTarget.getByteCount() +
            this.transparentColorTarget.getByteCount() +
            this.depthTargetTransparent.getByteCount() +
            (this.depthTargetOpaque
                ? this.depthTargetOpaque.getByteCount()
                : this.depthTextureOpaque.getByteCount()) +
            this.wboit.getByteCount() +
            this.dpoit.getByteCount() +
            this.marking.getByteCount() +
            this.postprocessing.getByteCount() +
            this.antialiasing.getByteCount() +
            this.dof.getByteCount()
        );
    }

    reset() {
        this.wboit.reset();
        this.dpoit.reset();
        this.postprocessing.reset();
        this.marking.invalidate();
    }

    setSize(width: number, height: number) {
        const w = this.colorTarget.getWidth();
        const h = this.colorTarget.getHeight();

        if (width !== w || height !== h) {
            this.colorTarget.setSize(width, height);
            this.depthTargetTransparent.setSize(width, height);
            this.transparentColorTarget.setSize(width, height);

            if (this.depthTargetOpaque) {
                this.depthTargetOpaque.setSize(width, height);
            } else {
                this.depthTextureOpaque.define(width, height);
            }

            ValueCell.update(this.copyFbo.values.uTexSize, Vec2.set(this.copyFbo.values.uTexSize.ref.value, width, height));
        }

        if (this.wboit.supported) {
            this.wboit.setSize(width, height);
        }

        if (this.dpoit.supported) {
            this.dpoit.setSize(width, height);
        }

        this.marking.setSize(width, height);
        this.postprocessing.setSize(width, height);
        this.antialiasing.setSize(width, height);
        this.dof.setSize(width, height);
    }

    renderEmissiveBloom(renderer: Renderer, camera: ICamera, scene: Scene, transparency: boolean): void {
        const bloom = this.postprocessing.bloom;
        bloom.emissiveTarget.bind();
        // (0,0,0,0) clear so the MAX blend in renderEmissiveTransparent isn't polluted by a white clear
        renderer.clear(false, true, true);
        // occlude emitters against real opaque depth so glow can't bleed through opaque foreground; packed depth builds its own
        const occludeWithOpaqueDepth = !this.packedDepth;
        if (occludeWithOpaqueDepth) {
            this.depthTextureOpaque.attachFramebuffer(bloom.emissiveTarget.framebuffer, 'depth');
        }
        renderer.renderEmissiveOpaque(scene.primitives, camera, this.depthTextureTransparent, occludeWithOpaqueDepth);
        if (transparency && scene.opacityAverage < 1) {
            renderer.renderEmissiveTransparent(scene.primitives, camera, this.depthTextureTransparent);
        }
        if (occludeWithOpaqueDepth) {
            this.depthTextureOpaque.detachFramebuffer(bloom.emissiveTarget.framebuffer, 'depth');
            bloom.emissiveTarget.depthRenderbuffer?.attachFramebuffer(bloom.emissiveTarget.framebuffer);
        }
    }

    private _renderBloom(renderer: Renderer, camera: ICamera, scene: Scene, postprocessingProps: PostprocessingProps): boolean {
        if (!BloomPass.isEnabled(postprocessingProps) || postprocessingProps.bloom.name !== 'on') return false;
        const bloom = this.postprocessing.bloom;
        const params = postprocessingProps.bloom.params;
        const emissiveBloom = params.mode === 'emissive';
        if (emissiveBloom && scene.emissiveAverage === 0) return false;

        if (emissiveBloom) {
            this.renderEmissiveBloom(renderer, camera, scene, params.transparency);
        }

        // Clear transparent buffers when opaque-only so luminosity doesn't sample stale halos.
        if (scene.opacityAverage >= 1) {
            this.transparentColorTarget.bind();
            renderer.clear(false, false, true);
            this.depthTargetTransparent.bind();
            renderer.clearDepth(true);
        }

        bloom.update(this.colorTarget.texture, this.transparentColorTarget.texture, bloom.emissiveTarget.texture, this.depthTextureOpaque, this.depthTextureTransparent, params, camera, true);
        bloom.render(camera.viewport);
        return true;
    }

    private _renderDpoit(renderer: Renderer, camera: ICamera, scene: Scene, iterations: number, transparentBackground: boolean, postprocessingProps: PostprocessingProps) {
        if (!this.dpoit.supported) throw new Error('expected dpoit to be supported');

        this.depthTextureOpaque.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
        renderer.clear(true);

        // render opaque primitives
        if (scene.hasOpaque) {
            renderer.renderOpaque(scene.primitives, camera);
        }

        this.depthTextureOpaque.detachFramebuffer(this.colorTarget.framebuffer, 'depth');

        if (PostprocessingPass.isTransparentDepthRequired(scene, postprocessingProps)) {
            this.depthTargetTransparent.bind();
            renderer.clearDepth(true);
            if (scene.opacityAverage < 1) {
                renderer.renderDepthTransparent(scene.primitives, camera, this.depthTextureOpaque);
            }
        }

        // render transparent primitives
        const isPostprocessingEnabled = PostprocessingPass.isEnabled(postprocessingProps);
        if (scene.opacityAverage < 1) {
            const target = isPostprocessingEnabled ? this.transparentColorTarget : this.colorTarget;
            if (isPostprocessingEnabled) {
                target.bind();
                renderer.clear(false, false, true);
            }

            const dpoitTextures = this.dpoit.bind();
            renderer.renderDpoitTransparent(scene.primitives, camera, this.depthTextureOpaque, dpoitTextures);

            for (let i = 0; i < iterations; i++) {
                if (isTimingMode) this.webgl.timer.mark('DpoitPass.layer');
                const dpoitTextures = this.dpoit.bindDualDepthPeeling();
                renderer.renderDpoitTransparent(scene.primitives, camera, this.depthTextureOpaque, dpoitTextures);

                if (iterations > 1) {
                    target.bind();
                    this.dpoit.renderBlendBack();
                }
                if (isTimingMode) this.webgl.timer.markEnd('DpoitPass.layer');
            }

            // evaluate dpoit
            target.bind();
            this.dpoit.render();
        }

        if (PostprocessingPass.isEnabled(postprocessingProps)) {
            this.postprocessing.render(camera, scene, false, transparentBackground, renderer.props.backgroundColor, postprocessingProps, renderer.light, renderer.ambientColor, this._renderBloom(renderer, camera, scene, postprocessingProps));
        }

        // render transparent volumes
        if (scene.volumes.renderables.length > 0) {
            renderer.renderVolume(scene.volumes, camera, this.depthTextureOpaque);
        }
    }

    private _renderWboit(renderer: Renderer, camera: ICamera, scene: Scene, transparentBackground: boolean, postprocessingProps: PostprocessingProps) {
        if (!this.wboit.supported) throw new Error('expected wboit to be supported');

        this.depthTextureOpaque.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
        renderer.clear(true);

        // render opaque primitives
        if (scene.hasOpaque) {
            renderer.renderOpaque(scene.primitives, camera);
        }

        if (PostprocessingPass.isTransparentDepthRequired(scene, postprocessingProps)) {
            this.depthTargetTransparent.bind();
            renderer.clearDepth(true);
            if (scene.opacityAverage < 1) {
                renderer.renderDepthTransparent(scene.primitives, camera, this.depthTextureOpaque);
            }
        }

        // render transparent primitives
        const isPostprocessingEnabled = PostprocessingPass.isEnabled(postprocessingProps);
        if (scene.opacityAverage < 1) {
            const target = isPostprocessingEnabled ? this.transparentColorTarget : this.colorTarget;
            if (isPostprocessingEnabled) {
                target.bind();
                renderer.clear(false, false, true);
            }

            this.wboit.bind();
            renderer.renderWboitTransparent(scene.primitives, camera, this.depthTextureOpaque);

            // evaluate wboit
            target.bind();
            this.wboit.render();
        }

        if (PostprocessingPass.isEnabled(postprocessingProps)) {
            this.postprocessing.render(camera, scene, false, transparentBackground, renderer.props.backgroundColor, postprocessingProps, renderer.light, renderer.ambientColor, this._renderBloom(renderer, camera, scene, postprocessingProps));
        }

        // render volumes
        if (scene.volumes.renderables.length > 0) {
            this.wboit.bind();
            renderer.renderWboitTransparent(scene.volumes, camera, this.depthTextureOpaque);

            // evaluate wboit
            const target = isPostprocessingEnabled ? this.postprocessing.target : this.colorTarget;
            target.bind();
            this.wboit.render();
        }

    }

    private _renderBlended(renderer: Renderer, camera: ICamera, scene: Scene, toDrawingBuffer: boolean, transparentBackground: boolean, postprocessingProps: PostprocessingProps) {
        if (toDrawingBuffer) {
            this.drawTarget.bind();
        } else {
            if (!this.packedDepth) {
                this.depthTextureOpaque.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
            } else {
                this.colorTarget.bind();
            }
        }

        renderer.clear(true);
        if (scene.hasOpaque) {
            renderer.renderOpaque(scene.primitives, camera);
        }

        if (!toDrawingBuffer) {
            // do a depth pass if not rendering to drawing buffer and
            // extensions.depthTexture is unsupported (i.e. depthTarget is set)
            if (this.depthTargetOpaque) {
                this.depthTargetOpaque.bind();
                renderer.clearDepth(true);
                renderer.renderDepthOpaque(scene.primitives, camera);
                this.colorTarget.bind();
            }

            if (PostprocessingPass.isTransparentDepthRequired(scene, postprocessingProps)) {
                this.depthTargetTransparent.bind();
                renderer.clearDepth(true);
                if (scene.opacityAverage < 1) {
                    renderer.renderDepthTransparent(scene.primitives, camera, this.depthTextureOpaque);
                }
            }

            // render transparent primitives
            const isPostprocessingEnabled = PostprocessingPass.isEnabled(postprocessingProps);
            if (scene.opacityAverage < 1) {
                if (isPostprocessingEnabled) {
                    this.transparentColorTarget.bind();
                    renderer.clear(false, false, true);

                    if (!this.packedDepth) {
                        this.depthTextureOpaque.attachFramebuffer(this.transparentColorTarget.framebuffer, 'depth');
                    } else {
                        this.colorTarget.depthRenderbuffer?.attachFramebuffer(this.transparentColorTarget.framebuffer);
                    }
                }

                renderer.renderBlendedTransparent(scene.primitives, camera);

                if (isPostprocessingEnabled) {
                    if (!this.packedDepth) {
                        this.depthTextureOpaque.detachFramebuffer(this.transparentColorTarget.framebuffer, 'depth');
                    } else {
                        this.colorTarget.depthRenderbuffer?.detachFramebuffer(this.transparentColorTarget.framebuffer);
                    }
                }
            }

            if (isPostprocessingEnabled) {
                if (!this.packedDepth) {
                    this.depthTextureOpaque.detachFramebuffer(this.postprocessing.target.framebuffer, 'depth');
                } else {
                    this.colorTarget.depthRenderbuffer?.detachFramebuffer(this.postprocessing.target.framebuffer);
                }

                this.postprocessing.render(camera, scene, false, transparentBackground, renderer.props.backgroundColor, postprocessingProps, renderer.light, renderer.ambientColor, this._renderBloom(renderer, camera, scene, postprocessingProps));

                if (!this.packedDepth) {
                    this.depthTextureOpaque.attachFramebuffer(this.postprocessing.target.framebuffer, 'depth');
                } else {
                    this.colorTarget.depthRenderbuffer?.attachFramebuffer(this.postprocessing.target.framebuffer);
                }
            }

            if (scene.volumes.renderables.length > 0) {
                const target = PostprocessingPass.isEnabled(postprocessingProps)
                    ? this.postprocessing.target : this.colorTarget;

                if (!this.packedDepth) {
                    this.depthTextureOpaque.detachFramebuffer(target.framebuffer, 'depth');
                } else {
                    this.colorTarget.depthRenderbuffer?.detachFramebuffer(target.framebuffer);
                }
                target.bind();

                renderer.renderVolume(scene.volumes, camera, this.depthTextureOpaque);

                if (!this.packedDepth) {
                    this.depthTextureOpaque.attachFramebuffer(target.framebuffer, 'depth');
                } else {
                    this.colorTarget.depthRenderbuffer?.attachFramebuffer(target.framebuffer);
                }
                target.bind();
            }
        } else if (scene.opacityAverage < 1) {
            renderer.renderBlendedTransparent(scene.primitives, camera);
        }
    }

    private _render(ctx: RenderContext<ICamera>, toDrawingBuffer: boolean, transparentBackground: boolean, props: Props, skipMarking: boolean) {
        const { renderer, camera, scene, frame } = ctx;
        if (camera.disabled) return;

        const volumeRendering = scene.volumes.renderables.length > 0;
        const postprocessingEnabled = PostprocessingPass.isEnabled(props.postprocessing);
        const antialiasingEnabled = AntialiasingPass.isEnabled(props.postprocessing);
        const dofEnabled = DofPass.isEnabled(props.postprocessing);
        const hasMarking = !skipMarking && MarkingPass.hasMarking(scene, props);
        // marking is composed onto the finished image, which therefore has to be readable
        const offscreen = !toDrawingBuffer || hasMarking;

        const { x, y, width, height } = camera.viewport;
        renderer.setViewport(x, y, width, height);
        renderer.update(camera, scene, frame);

        if (transparentBackground && !antialiasingEnabled && !dofEnabled && !offscreen && !postprocessingEnabled) {
            this.drawTarget.bind();
            renderer.clear(false);
        }

        let oitEnabled = false;
        if (this.transparencyMode === 'wboit' && this.wboit.supported) {
            this._renderWboit(renderer, camera, scene, transparentBackground, props.postprocessing);
            oitEnabled = true;
        } else if (this.transparencyMode === 'dpoit' && this.dpoit.supported) {
            this._renderDpoit(renderer, camera, scene, props.dpoitIterations, transparentBackground, props.postprocessing);
            oitEnabled = true;
        } else {
            this._renderBlended(renderer, camera, scene, !volumeRendering && !postprocessingEnabled && !antialiasingEnabled && !dofEnabled && !offscreen, transparentBackground, props.postprocessing);
        }

        const target = postprocessingEnabled
            ? this.postprocessing.target
            : offscreen || volumeRendering || oitEnabled || dofEnabled
                ? this.colorTarget
                : this.drawTarget;

        target.bind();
        this._renderHelpers(ctx, target);

        const output = this._renderOutput(camera, scene, props, { offscreen, postprocessingEnabled, antialiasingEnabled, dofEnabled, volumeRendering, oitEnabled });

        if (!output) {
            this.marking.invalidate();
        } else if (!skipMarking) {
            this.marking.present(ctx, props, { base: output, toDrawingBuffer, offsets: SingleSample, restart: true, samples: 1, shading: this.getMarkingShading(props.postprocessing) });
        }

        this.webgl.gl.flush();
    }

    private _renderHelpers(ctx: RenderContext<ICamera>, target: RenderTarget) {
        const { renderer, camera, helper, frame } = ctx;
        if (helper.debug.isEnabled || helper.pointer.isEnabled) {
            if (!this.packedDepth) {
                this.depthTextureOpaque.attachFramebuffer(target.framebuffer, 'depth');
            }
            if (helper.debug.isEnabled) {
                helper.debug.syncVisibility();
                for (const scene of helper.debug.scenes) {
                    renderer.renderBlended(scene, camera);
                }
            }
            if (helper.pointer.isEnabled) {
                helper.pointer.setCamera(camera);
                renderer.update(helper.pointer.camera, helper.pointer.scene, frame);
                renderer.renderBlended(helper.pointer.scene, helper.pointer.camera);
            }
            if (!this.packedDepth) {
                this.depthTextureOpaque.detachFramebuffer(target.framebuffer, 'depth');
            }
        }
        if (helper.handle.isEnabled) {
            renderer.renderBlended(helper.handle.scene, camera);
        }
        if (helper.camera.isEnabled) {
            helper.camera.update(camera);
            renderer.update(helper.camera.camera, helper.camera.scene, frame);
            renderer.renderBlended(helper.camera.scene, helper.camera.camera);
        }
    }

    /** Returns the target the finished image ended up in, or undefined when it went to the drawing buffer. */
    private _renderOutput(camera: ICamera, scene: Scene, props: Props, flags: { offscreen: boolean, postprocessingEnabled: boolean, antialiasingEnabled: boolean, dofEnabled: boolean, volumeRendering: boolean, oitEnabled: boolean }): RenderTarget | undefined {
            const { offscreen, postprocessingEnabled, antialiasingEnabled, dofEnabled, volumeRendering, oitEnabled } = flags;
            const toBuffer = !offscreen;

        let needsTargetCopy = false;

        if (antialiasingEnabled) {
            const input = postprocessingEnabled
                ? this.postprocessing.target.texture
                : this.colorTarget.texture;
            this.antialiasing.render(camera, input, toBuffer && !dofEnabled, props.postprocessing);
        } else if (toBuffer && !dofEnabled) {
            needsTargetCopy = true;
        }

        if (dofEnabled && props.postprocessing.dof.name === 'on') {
            const input = antialiasingEnabled
                ? this.antialiasing.target.texture
                : postprocessingEnabled
                    ? this.postprocessing.target.texture
                    : this.colorTarget.texture;
            this.dof.update(camera, input, this.depthTargetOpaque?.texture || this.depthTextureOpaque, this.depthTextureTransparent, props.postprocessing.dof.params, scene.boundingSphereVisible);
            this.dof.render(camera.viewport, toBuffer ? undefined : this.getColorTarget(props.postprocessing));
        } else if (toBuffer && !antialiasingEnabled) {
            needsTargetCopy = true;
        }

        if (needsTargetCopy) {
            if (postprocessingEnabled) {
                this.copyToDrawingBuffer(this.postprocessing.target, camera.viewport);
            } else if (volumeRendering || oitEnabled) {
                this.copyToDrawingBuffer(this.colorTarget, camera.viewport);
            }
        }

        return offscreen ? this.getColorTarget(props.postprocessing) : undefined;
    }

    private copyToDrawingBuffer(src: RenderTarget, viewport: Viewport) {
        const { gl, state } = this.webgl;

        if (this.copyFbo.values.tColor.ref.value !== src.texture) {
            ValueCell.update(this.copyFbo.values.tColor, src.texture);
            this.copyFbo.update();
        }

        this.drawTarget.bind();
        state.enable(gl.SCISSOR_TEST);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);

        const { x, y, width, height } = viewport;
        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);
        this.copyFbo.render();
    }

    render(ctx: RenderContext, props: Props, toDrawingBuffer: boolean, skipMarking = false) {
        if (isTimingMode) this.webgl.timer.mark('DrawPass.render');
        const { renderer, camera, scene } = ctx;

        this.postprocessing.setTransparentBackground(props.transparentBackground);
        const transparentBackground = this.isTransparentBackground(scene, props);

        renderer.setTransparentBackground(transparentBackground);
        renderer.setDrawingBufferSize(this.colorTarget.getWidth(), this.colorTarget.getHeight());
        renderer.setPixelRatio(this.webgl.pixelRatio);

        if (StereoCamera.is(camera)) {
            if (isTimingMode) this.webgl.timer.mark('StereoCamera.left');
            this._render({ ...ctx, camera: camera.left }, toDrawingBuffer, transparentBackground, props, skipMarking);
            if (isTimingMode) this.webgl.timer.markEnd('StereoCamera.left');
            if (isTimingMode) this.webgl.timer.mark('StereoCamera.right');
            this._render({ ...ctx, camera: camera.right }, toDrawingBuffer, transparentBackground, props, skipMarking);
            if (isTimingMode) this.webgl.timer.markEnd('StereoCamera.right');
            // only the second eye is left in the target, which is not a valid recompose source
            this.marking.invalidate();
        } else {
            this._render({ ...ctx, camera }, toDrawingBuffer, transparentBackground, props, skipMarking);
        }
        if (isTimingMode) this.webgl.timer.markEnd('DrawPass.render');
    }

    private isTransparentBackground(scene: Scene, props: Props) {
        const pp = props.postprocessing;
        const backgroundEnabled = this.postprocessing.background.isEnabled(pp);
        // premultiplied scene so bloom can composite the solid background last
        const bloomCompositesBackground = !props.transparentBackground && !backgroundEnabled &&
            BloomPass.isEnabled(pp) && pp.bloom.name === 'on' &&
            !(pp.bloom.params.mode === 'emissive' && scene.emissiveAverage === 0);
        return props.transparentBackground || backgroundEnabled || bloomCompositesBackground;
    }

    /** shading of the last render that marking tint and dim keep, or null if there is none */
    getMarkingShading(postprocessingProps: PostprocessingProps): MarkingShading | null {
        return SsaoPass.isEnabled(postprocessingProps) ? this.ssaoShading : null;
    }

    getColorTarget(postprocessingProps: PostprocessingProps): RenderTarget {
        if (DofPass.isEnabled(postprocessingProps)) {
            return this.dof.target;
        } else if (AntialiasingPass.isEnabled(postprocessingProps)) {
            return this.antialiasing.target;
        } else if (PostprocessingPass.isEnabled(postprocessingProps)) {
            return this.postprocessing.target;
        }
        return this.colorTarget;
    }
}
