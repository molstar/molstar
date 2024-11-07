/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { Renderer, RendererParams } from '../../mol-gl/renderer';
import { Scene } from '../../mol-gl/scene';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { DrawPass } from './draw';
import { PostprocessingParams } from './postprocessing';
import { MultiSamplePass, MultiSampleParams, MultiSampleHelper } from './multi-sample';
import { Camera } from '../camera';
import { Viewport } from '../camera/util';
import { PixelData } from '../../mol-util/image';
import { Helper } from '../helper/helper';
import { CameraHelper, CameraHelperParams } from '../helper/camera-helper';
import { MarkingParams } from './marking';
import { AssetManager } from '../../mol-util/assets';
import { IlluminationParams, IlluminationPass } from './illumination';
import { RuntimeContext } from '../../mol-task';
import { isTimingMode } from '../../mol-util/debug';
import { printTimerResults } from '../../mol-gl/webgl/timer';

export const ImageParams = {
    transparentBackground: PD.Boolean(false),
    dpoitIterations: PD.Numeric(2, { min: 1, max: 10, step: 1 }),
    multiSample: PD.Group(MultiSampleParams),
    postprocessing: PD.Group(PostprocessingParams),
    marking: PD.Group(MarkingParams),
    illumination: PD.Group(IlluminationParams),

    cameraHelper: PD.Group(CameraHelperParams),
    renderer: PD.Group(RendererParams),
};
export type ImageProps = PD.Values<typeof ImageParams>

export class ImagePass {
    private _width = 0;
    private _height = 0;
    private _camera = new Camera();

    readonly props: ImageProps;

    private _colorTarget: RenderTarget;
    get colorTarget() { return this._colorTarget; }

    private readonly drawPass: DrawPass;
    private readonly illuminationPass: IlluminationPass;
    private readonly multiSamplePass: MultiSamplePass;
    private readonly multiSampleHelper: MultiSampleHelper;
    private readonly helper: Helper;

    get width() { return this._width; }
    get height() { return this._height; }

    constructor(private webgl: WebGLContext, assetManager: AssetManager, private renderer: Renderer, private scene: Scene, private camera: Camera, helper: Helper, transparency: 'wboit' | 'dpoit' | 'blended', props: Partial<ImageProps>) {
        this.props = { ...PD.getDefaultValues(ImageParams), ...props };

        this.drawPass = new DrawPass(webgl, assetManager, 128, 128, transparency);
        this.illuminationPass = new IlluminationPass(webgl, this.drawPass);
        this.multiSamplePass = new MultiSamplePass(webgl, this.drawPass);
        this.multiSampleHelper = new MultiSampleHelper(this.multiSamplePass);

        this.helper = {
            camera: new CameraHelper(webgl, this.props.cameraHelper),
            debug: helper.debug,
            handle: helper.handle,
        };

        this.setSize(1024, 768);
    }

    updateBackground() {
        return new Promise<void>(resolve => {
            this.drawPass.postprocessing.background.update(this.camera, this.props.postprocessing.background, () => {
                resolve();
            });
        });
    }

    setSize(width: number, height: number) {
        if (width === this._width && height === this._height) return;

        this._width = width;
        this._height = height;

        this.drawPass.setSize(width, height);
        this.illuminationPass.setSize(width, height);
        this.multiSamplePass.syncSize();
    }

    setProps(props: Partial<ImageProps> = {}) {
        Object.assign(this.props, props);
        if (props.cameraHelper) this.helper.camera.setProps(props.cameraHelper);
    }

    async render(runtime: RuntimeContext) {
        Camera.copySnapshot(this._camera.state, this.camera.state);
        Viewport.set(this._camera.viewport, 0, 0, this._width, this._height);
        this._camera.update();

        const ctx = { renderer: this.renderer, camera: this._camera, scene: this.scene, helper: this.helper };
        if (this.illuminationPass.supported && this.props.illumination.enabled) {
            await runtime.update({ message: 'Tracing...', current: 1, max: this.illuminationPass.getMaxIterations(this.props) });
            this.illuminationPass.reset(true);
            while (this.illuminationPass.shouldRender(this.props)) {
                if (isTimingMode) this.webgl.timer.mark('ImagePass.render', { captureStats: true });
                this.illuminationPass.render(ctx, this.props, false);
                if (isTimingMode) this.webgl.timer.markEnd('ImagePass.render');
                if (runtime.shouldUpdate) {
                    await runtime.update({ current: this.illuminationPass.iteration });
                }
                await this.webgl.waitForGpuCommandsComplete();
            }
            this._colorTarget = this.illuminationPass.colorTarget;
        } else {
            if (isTimingMode) this.webgl.timer.mark('ImagePass.render', { captureStats: true });
            if (MultiSamplePass.isEnabled(this.props.multiSample)) {
                this.multiSampleHelper.render(ctx, this.props, false);
                this._colorTarget = this.multiSamplePass.colorTarget;
            } else {
                this.drawPass.render(ctx, this.props, false);
                this._colorTarget = this.drawPass.getColorTarget(this.props.postprocessing);
            }
            if (isTimingMode) this.webgl.timer.markEnd('ImagePass.render');
        }

        if (isTimingMode) {
            const timerResults = this.webgl.timer.resolve();
            if (timerResults) {
                for (const result of timerResults) {
                    printTimerResults([result]);
                }
            }
        }

        if (isTimingMode) {
            const timerResults = this.webgl.timer.resolve();
            if (timerResults) {
                for (const result of timerResults) {
                    printTimerResults([result]);
                }
            }
        }
    }

    async getImageData(runtime: RuntimeContext, width: number, height: number, viewport?: Viewport) {
        this.setSize(width, height);
        await this.render(runtime);
        this.colorTarget.bind();

        const w = viewport?.width ?? width, h = viewport?.height ?? height;

        const array = new Uint8Array(w * h * 4);
        if (!viewport) {
            this.webgl.readPixels(0, 0, w, h, array);
        } else {
            this.webgl.readPixels(viewport.x, height - viewport.y - viewport.height, w, h, array);
        }
        const pixelData = PixelData.create(array, w, h);
        PixelData.flipY(pixelData);
        PixelData.divideByAlpha(pixelData);
        return new ImageData(new Uint8ClampedArray(array), w, h);
    }
}