import * as HME from 'h264-mp4-encoder';
import { Viewport } from '../../mol-canvas3d/camera/util';
import { AnimateCameraSpin } from '../../mol-plugin-state/animation/built-in/camera-spin';
import { PluginContext } from '../../mol-plugin/context';

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

    sleep() {
        return new Promise(res => setTimeout(res, 16.6));
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
        await this.plugin.managers.animation.play(AnimateCameraSpin, { durationInMs: durationMs, direction: 'cw' });

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