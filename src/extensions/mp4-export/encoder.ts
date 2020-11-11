import * as HME from 'h264-mp4-encoder';
import { Viewport } from '../../mol-canvas3d/camera/util';
import { ImagePass } from '../../mol-canvas3d/passes/image';
import { AnimateCameraSpin } from '../../mol-plugin-state/animation/built-in/camera-spin';
import { PluginStateAnimation } from '../../mol-plugin-state/animation/model';
import { PluginContext } from '../../mol-plugin/context';
import { RuntimeContext } from '../../mol-task';
import { Color } from '../../mol-util/color';

export interface Mp4EncoderParams<A extends PluginStateAnimation = PluginStateAnimation> {
    pass: ImagePass,
    customBackground?: Color,
    animation: PluginStateAnimation.Instance<A>,
    width: number,
    height: number,
    viewport?: Viewport,
    /** default is 30 */
    fps?: number,
    /** Number from 10 (best quality, slowest) to 51 (worst, fastest) */
    quantizationParameter?: number
}

export async function encodeMp4Animation<A extends PluginStateAnimation>(plugin: PluginContext, ctx: RuntimeContext, params: Mp4EncoderParams<A>) {
    await ctx.update({ message: 'Initializing...', isIndeterminate: true });

    validateViewport(params);
    const durationMs = PluginStateAnimation.getDuration(plugin, params.animation);
    if (durationMs === void 0) {
        throw new Error('The animation does not have the duration specified.');
    }

    const encoder = await HME.createH264MP4Encoder();

    const { width, height } = params;
    const vw = params.viewport?.width ?? width, vh = params.viewport?.height ?? height;

    encoder.width = vw;
    encoder.height = vh;
    if (params.quantizationParameter) encoder.quantizationParameter = params.quantizationParameter;
    if (params.fps) encoder.frameRate = params.fps;
    encoder.initialize();

    const loop = plugin.animationLoop;
    const originalBackground = params.customBackground ? plugin.canvas3d?.props.renderer.backgroundColor : void 0;
    let stoppedAnimation = true, finalized = false;

    try {
        loop.stop();
        loop.resetTime(0);
        plugin.canvas3d?.setProps({ renderer: { backgroundColor: params.customBackground } }, true);

        const fps = encoder.frameRate;
        const N = Math.ceil(durationMs / 1000 * fps);
        const dt = durationMs / N;

        await ctx.update({ message: 'Rendering...', isIndeterminate: false, current: 0, max: N + 1 });

        await plugin.managers.animation.play(params.animation.definition, params.animation.params);
        stoppedAnimation = false;
        for (let i = 0; i <= N; i++) {
            await loop.tick(i * dt, { isSynchronous: true, manualDraw: true });

            const image = params.pass.getImageData(width, height, params.viewport);
            encoder.addFrameRgba(image.data);

            if (ctx.shouldUpdate) {
                await ctx.update({ current: i + 1 });
            }
        }
        await ctx.update({ message: 'Applying finishing touches...', isIndeterminate: true });
        await plugin.managers.animation.stop();
        stoppedAnimation = true;
        encoder.finalize();
        finalized = true;
        return encoder.FS.readFile(encoder.outputFilename);
    } finally {
        if (finalized) encoder.delete();
        if (originalBackground) {
            plugin.canvas3d?.setProps({ renderer: { backgroundColor: originalBackground } }, true);
        }
        if (!stoppedAnimation) await plugin.managers.animation.stop();
        loop.start();
    }
}

function validateViewport(params: Mp4EncoderParams) {
    if (!params.viewport) return;

    if (params.viewport.x + params.viewport.width > params.width || params.viewport.x + params.viewport.width >= params.width) {
        throw new Error('Viewport exceeds the canvas dimensions.');
    }
}


export class Mp4Encoder {
    createImagePass() {
        const pass = this.plugin.canvas3d!.getImagePass({
            transparentBackground: true,
            cameraHelper: { axes: { name: 'off', params: {} } },
            multiSample: { mode: 'on', sampleLevel: 4 }, // { mode: 'on', sampleLevel: 2 },
            postprocessing: this.plugin.canvas3d!.props.postprocessing
        });
        return pass;
    }

    async generate() {
        // const w = 1024, h = 768;
        // const w = 1920, h = 1080;
        const w = 1280, h = 720;

        const viewport: Viewport = {
            x: w / 2 - 300,
            y: h / 2 - 250,
            width: 600,
            height: 500
        };

        const encoder = await HME.createH264MP4Encoder();
        encoder.width = viewport.width;
        encoder.height = viewport.height;
        encoder.frameRate = 30;
        encoder.quantizationParameter = 15;
        encoder.initialize();

        console.log('creating image pass');
        const pass = this.createImagePass();

        // const canvas = document.createElement('canvas');
        // canvas.width = w;
        // canvas.height = h;
        // const canvasCtx = canvas.getContext('2d')!;

        const loop = this.plugin.animationLoop;

        loop.stop();
        loop.resetTime(0);

        const durationMs = 2500;
        const fps = encoder.frameRate;

        const color = this.plugin.canvas3d?.props.renderer.backgroundColor;
        this.plugin.canvas3d?.setProps({ renderer: { backgroundColor: 0x000000 as any } }, true);

        // await this.plugin.managers.animation.play(AnimateAssemblyUnwind, { durationInMs: 3000, playOnce: true, target: 'all' });
        await this.plugin.managers.animation.play(AnimateCameraSpin, { durationInMs: durationMs, speed: 1, direction: 'cw' });

        // const imageData: Uint8ClampedArray[] = [];
        const N = Math.ceil(durationMs / 1000 * fps);
        const dt = durationMs / N;
        for (let i = 0; i <= N; i++) {
            await loop.tick(i * dt, { isSynchronous: true, manualDraw: true });

            const image = pass.getImageData(w, h, viewport);
            encoder.addFrameRgba(image.data);
            // if (i === 0) canvasCtx.putImageData(image, 0, 0);

            console.log(`frame ${i + 1}/${N + 1}`);
        }

        this.plugin.canvas3d?.setProps({ renderer: { backgroundColor: color } }, true);

        console.log('finalizing');
        encoder.finalize();
        console.log('finalized');
        const uint8Array = encoder.FS.readFile(encoder.outputFilename);
        console.log('encoded');
        encoder.delete();

        await this.plugin.managers.animation.stop();
        loop.start();

        return { movie: uint8Array };
    }

    constructor(private plugin: PluginContext) {
    }
}