/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { WebGPUContext } from '../../mol-gl/webgpu/context';
import { copyWebGPUCameraState } from '../../mol-gl/webgpu/camera';
import { WebGPURenderer } from '../../mol-gl/webgpu/renderer';
import { AssetManager } from '../../mol-util/assets';
import { RuntimeContext } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Camera } from '../camera';
import { Viewport } from '../camera/util';
import { ImageParams, ImageProps } from './image';
import { createWebGPUCameraHelper } from '../helper/camera-helper';
import { createWebGPUHandleHelper } from '../helper/handle-helper';
import { createWebGPUPointerHelper } from '../helper/pointer-helper';
import { WebGPUDebugRegistry } from '../helper/debug-registry';

export class WebGPUImagePass {
    readonly props: ImageProps;
    private readonly renderer: Promise<WebGPURenderer>;
    private readonly camera = new Camera();
    private readonly cameraHelper: ReturnType<typeof createWebGPUCameraHelper>;
    private _width = 1024;
    private _height = 768;
    get width() { return this._width; }
    get height() { return this._height; }

    constructor(private readonly context: WebGPUContext, private readonly sourceCamera: Camera, private readonly objects: () => GraphicsRenderObject[], props: Partial<ImageProps>, assets?: AssetManager,
        private readonly worldHelpers?: { handle?: ReturnType<typeof createWebGPUHandleHelper>, pointer?: ReturnType<typeof createWebGPUPointerHelper>, debug?: WebGPUDebugRegistry, transparency?: () => 'blended' | 'wboit' | 'dpoit' }) {
        this.props = PD.merge(ImageParams, PD.getDefaultValues(ImageParams), props);
        this.cameraHelper = createWebGPUCameraHelper(() => 1, this.props.cameraHelper);
        this.renderer = WebGPURenderer.create(context, assets);
    }

    getByteCount() { return this.width * this.height * 24; }
    async updateBackground() { await (await this.renderer).updateBackground(this.props.postprocessing); }
    setSize(width: number, height: number) { this._width = width; this._height = height; }
    setProps(props: Partial<ImageProps> = {}) { Object.assign(this.props, PD.merge(ImageParams, this.props, props)); if (props.cameraHelper) this.cameraHelper.setProps(this.props.cameraHelper); }

    async render(runtime: RuntimeContext) {
        copyWebGPUCameraState(this.camera, this.sourceCamera);
        Object.assign(this.camera.viewport, { x: 0, y: 0, width: this.width, height: this.height });
        this.camera.update();
        const renderer = await this.renderer;
        renderer.setTransparency(this.worldHelpers?.transparency?.() ?? 'blended');
        renderer.setDpoitIterations(this.props.dpoitIterations);
        renderer.setCameraHelper(this.cameraHelper);
        renderer.setWorldHelpers(this.worldHelpers?.handle, this.worldHelpers?.pointer);
        renderer.setDebugRegistry(this.worldHelpers?.debug);
        await renderer.updateBackground(this.props.postprocessing);
        if (this.props.illumination.enabled) {
            const max = Math.pow(2, Math.max(0, Math.min(16, Math.round(this.props.illumination.maxIterations))));
            for (let iteration = 0; iteration < max; iteration++) {
                if (runtime.shouldUpdate) await runtime.update({ message: 'Tracing...', current: iteration + 1, max });
                renderer.render(this.objects(), this.camera, this.props.renderer, this.props.transparentBackground, 1, { width: this.width, height: this.height, present: false }, this.props.postprocessing, this.props.marking, this.props.multiSample, iteration === 0, false, this.props.illumination);
                await this.context.device.queue.onSubmittedWorkDone();
            }
        } else {
            renderer.render(this.objects(), this.camera, this.props.renderer, this.props.transparentBackground, 1, { width: this.width, height: this.height, present: false }, this.props.postprocessing, this.props.marking, this.props.multiSample, true, true);
        }
    }

    async getImageRaw(runtime: RuntimeContext, width: number, height: number, viewport?: Viewport) {
        this.setSize(width, height);
        await this.render(runtime);
        const pixels = await (await this.renderer).readPixels();
        const w = viewport?.width ?? width, h = viewport?.height ?? height;
        const x = viewport?.x ?? 0, y = viewport?.y ?? 0;
        const array = new Uint8ClampedArray(w * h * 4);
        for (let row = 0; row < h; row++) {
            array.set(pixels.array.subarray(((y + row) * width + x) * 4, ((y + row) * width + x + w) * 4), row * w * 4);
        }
        // GPU colors are premultiplied, whereas ImageData uses straight alpha.
        for (let i = 0; i < array.length; i += 4) {
            const alpha = array[i + 3];
            if (alpha > 0 && alpha < 255) for (let c = 0; c < 3; c++) array[i + c] = array[i + c] * 255 / alpha;
        }
        return { data: array, width: w, height: h };
    }

    async getImageData(runtime: RuntimeContext, width: number, height: number, viewport?: Viewport) {
        const image = await this.getImageRaw(runtime, width, height, viewport);
        return new ImageData(image.data, image.width, image.height);
    }

    async dispose() { (await this.renderer).dispose(); this.cameraHelper.scene.clear(); }
}
