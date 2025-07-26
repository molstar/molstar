/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PickingId } from '../../mol-geo/geometry/picking';
import { PickType, Renderer } from '../../mol-gl/renderer';
import { Scene } from '../../mol-gl/scene';
import { PixelPackBuffer } from '../../mol-gl/webgl/buffer';
import { isWebGL2 } from '../../mol-gl/webgl/compat';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Framebuffer } from '../../mol-gl/webgl/framebuffer';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { Renderbuffer } from '../../mol-gl/webgl/renderbuffer';
import { Texture } from '../../mol-gl/webgl/texture';
import { Vec3 } from '../../mol-math/linear-algebra';
import { isDebugMode, isTimingMode } from '../../mol-util/debug';
import { now } from '../../mol-util/now';
import { unpackRGBAToDepth, unpackRGBToInt } from '../../mol-util/number-packing';
import { ICamera } from '../camera';
import { Viewport } from '../camera/util';
import { Helper } from '../helper/helper';

const NullId = Math.pow(2, 24) - 2;

export type PickData = { id: PickingId, position: Vec3 }

export type AsyncPickData = {
    tryGet: () => 'pending' | PickData | undefined,
}

export const DefaultPickOptions = {
    pickPadding: 1,
    maxAsyncReadLag: 5,
};
export type PickOptions = typeof DefaultPickOptions

//

export class PickPass {
    private readonly objectPickTarget: RenderTarget;
    private readonly instancePickTarget: RenderTarget;
    private readonly groupPickTarget: RenderTarget;
    private readonly depthPickTarget: RenderTarget;

    private readonly framebuffer: Framebuffer;

    private readonly objectPickTexture: Texture;
    private readonly instancePickTexture: Texture;
    private readonly groupPickTexture: Texture;
    private readonly depthPickTexture: Texture;

    private readonly objectPickFramebuffer: Framebuffer;
    private readonly instancePickFramebuffer: Framebuffer;
    private readonly groupPickFramebuffer: Framebuffer;
    private readonly depthPickFramebuffer: Framebuffer;

    private readonly depthRenderbuffer: Renderbuffer;

    private pickWidth: number;
    private pickHeight: number;

    constructor(private webgl: WebGLContext, private width: number, private height: number, private pickScale: number) {
        const pickRatio = pickScale / webgl.pixelRatio;
        this.pickWidth = Math.ceil(width * pickRatio);
        this.pickHeight = Math.ceil(height * pickRatio);

        const { resources, extensions: { drawBuffers }, gl } = webgl;

        if (drawBuffers) {
            this.objectPickTexture = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
            this.objectPickTexture.define(this.pickWidth, this.pickHeight);

            this.instancePickTexture = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
            this.instancePickTexture.define(this.pickWidth, this.pickHeight);

            this.groupPickTexture = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
            this.groupPickTexture.define(this.pickWidth, this.pickHeight);

            this.depthPickTexture = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
            this.depthPickTexture.define(this.pickWidth, this.pickHeight);

            this.framebuffer = resources.framebuffer();

            this.objectPickFramebuffer = resources.framebuffer();
            this.instancePickFramebuffer = resources.framebuffer();
            this.groupPickFramebuffer = resources.framebuffer();
            this.depthPickFramebuffer = resources.framebuffer();

            this.framebuffer.bind();
            drawBuffers!.drawBuffers([
                drawBuffers!.COLOR_ATTACHMENT0,
                drawBuffers!.COLOR_ATTACHMENT1,
                drawBuffers!.COLOR_ATTACHMENT2,
                drawBuffers!.COLOR_ATTACHMENT3,
            ]);

            this.objectPickTexture.attachFramebuffer(this.framebuffer, 'color0');
            this.instancePickTexture.attachFramebuffer(this.framebuffer, 'color1');
            this.groupPickTexture.attachFramebuffer(this.framebuffer, 'color2');
            this.depthPickTexture.attachFramebuffer(this.framebuffer, 'color3');

            this.depthRenderbuffer = isWebGL2(gl)
                ? resources.renderbuffer('depth32f', 'depth', this.pickWidth, this.pickHeight)
                : resources.renderbuffer('depth16', 'depth', this.pickWidth, this.pickHeight);

            this.depthRenderbuffer.attachFramebuffer(this.framebuffer);

            this.objectPickTexture.attachFramebuffer(this.objectPickFramebuffer, 'color0');
            this.instancePickTexture.attachFramebuffer(this.instancePickFramebuffer, 'color0');
            this.groupPickTexture.attachFramebuffer(this.groupPickFramebuffer, 'color0');
            this.depthPickTexture.attachFramebuffer(this.depthPickFramebuffer, 'color0');
        } else {
            this.objectPickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
            this.instancePickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
            this.groupPickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
            this.depthPickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
        }
    }

    dispose() {
        if (this.webgl.extensions.drawBuffers) {
            this.framebuffer.destroy();

            this.objectPickTexture.destroy();
            this.instancePickTexture.destroy();
            this.groupPickTexture.destroy();
            this.depthPickTexture.destroy();

            this.objectPickFramebuffer.destroy();
            this.instancePickFramebuffer.destroy();
            this.groupPickFramebuffer.destroy();
            this.depthPickFramebuffer.destroy();

            this.depthRenderbuffer.destroy();
        } else {
            this.objectPickTarget.destroy();
            this.instancePickTarget.destroy();
            this.groupPickTarget.destroy();
            this.depthPickTarget.destroy();
        }
    }

    get pickRatio() {
        return this.pickScale / this.webgl.pixelRatio;
    }

    setPickScale(pickScale: number) {
        this.pickScale = pickScale;
        this.setSize(this.width, this.height);
    }

    bindObject() {
        if (this.webgl.extensions.drawBuffers) {
            this.objectPickFramebuffer.bind();
        } else {
            this.objectPickTarget.bind();
        }
    }

    bindInstance() {
        if (this.webgl.extensions.drawBuffers) {
            this.instancePickFramebuffer.bind();
        } else {
            this.instancePickTarget.bind();
        }
    }

    bindGroup() {
        if (this.webgl.extensions.drawBuffers) {
            this.groupPickFramebuffer.bind();
        } else {
            this.groupPickTarget.bind();
        }
    }

    bindDepth() {
        if (this.webgl.extensions.drawBuffers) {
            this.depthPickFramebuffer.bind();
        } else {
            this.depthPickTarget.bind();
        }
    }

    get drawingBufferHeight() {
        return this.height;
    }

    setSize(width: number, height: number) {
        this.width = width;
        this.height = height;

        const pickRatio = this.pickScale / this.webgl.pixelRatio;
        const pickWidth = Math.ceil(this.width * pickRatio);
        const pickHeight = Math.ceil(this.height * pickRatio);

        if (pickWidth !== this.pickWidth || pickHeight !== this.pickHeight) {
            this.pickWidth = pickWidth;
            this.pickHeight = pickHeight;

            if (this.webgl.extensions.drawBuffers) {
                this.objectPickTexture.define(this.pickWidth, this.pickHeight);
                this.instancePickTexture.define(this.pickWidth, this.pickHeight);
                this.groupPickTexture.define(this.pickWidth, this.pickHeight);
                this.depthPickTexture.define(this.pickWidth, this.pickHeight);

                this.depthRenderbuffer.setSize(this.pickWidth, this.pickHeight);
            } else {
                this.objectPickTarget.setSize(this.pickWidth, this.pickHeight);
                this.instancePickTarget.setSize(this.pickWidth, this.pickHeight);
                this.groupPickTarget.setSize(this.pickWidth, this.pickHeight);
                this.depthPickTarget.setSize(this.pickWidth, this.pickHeight);
            }
        }
    }

    reset() {
        const { drawBuffers } = this.webgl.extensions;

        if (drawBuffers) {
            this.framebuffer.bind();
            drawBuffers!.drawBuffers([
                drawBuffers!.COLOR_ATTACHMENT0,
                drawBuffers!.COLOR_ATTACHMENT1,
                drawBuffers!.COLOR_ATTACHMENT2,
                drawBuffers!.COLOR_ATTACHMENT3,
            ]);

            this.objectPickTexture.attachFramebuffer(this.framebuffer, 'color0');
            this.instancePickTexture.attachFramebuffer(this.framebuffer, 'color1');
            this.groupPickTexture.attachFramebuffer(this.framebuffer, 'color2');
            this.depthPickTexture.attachFramebuffer(this.framebuffer, 'color3');

            this.depthRenderbuffer.attachFramebuffer(this.framebuffer);

            this.objectPickTexture.attachFramebuffer(this.objectPickFramebuffer, 'color0');
            this.instancePickTexture.attachFramebuffer(this.instancePickFramebuffer, 'color0');
            this.groupPickTexture.attachFramebuffer(this.groupPickFramebuffer, 'color0');
            this.depthPickTexture.attachFramebuffer(this.depthPickFramebuffer, 'color0');
        }
    }

    private renderVariant(renderer: Renderer, camera: ICamera, scene: Scene, helper: Helper, variant: 'pick' | 'depth', pickType: number) {
        renderer.clear(false);
        renderer.update(camera, scene);
        renderer.renderPick(scene.primitives, camera, variant, pickType);

        if (helper.handle.isEnabled) {
            renderer.renderPick(helper.handle.scene, camera, variant, pickType);
        }

        if (helper.camera.isEnabled) {
            helper.camera.update(camera);
            renderer.update(helper.camera.camera, helper.camera.scene);
            renderer.renderPick(helper.camera.scene, helper.camera.camera, variant, pickType);
        }
    }

    render(renderer: Renderer, camera: ICamera, scene: Scene, helper: Helper) {
        if (this.webgl.extensions.drawBuffers) {
            this.framebuffer.bind();
            this.renderVariant(renderer, camera, scene, helper, 'pick', PickType.None);
            // printTextureImage(readTexture(this.webgl, this.groupPickTexture, new Uint8Array(this.pickWidth * this.pickHeight * 4)), { scale: 16, id: 'group', pixelated: true, useCanvas: true, flipY: true });
        } else {
            this.objectPickTarget.bind();
            this.renderVariant(renderer, camera, scene, helper, 'pick', PickType.Object);

            this.instancePickTarget.bind();
            this.renderVariant(renderer, camera, scene, helper, 'pick', PickType.Instance);

            this.groupPickTarget.bind();
            this.renderVariant(renderer, camera, scene, helper, 'pick', PickType.Group);
            // printTextureImage(readTexture(this.webgl, this.groupPickTarget.texture, new Uint8Array(this.pickWidth * this.pickHeight * 4)), { scale: 16, id: 'group', pixelated: true, useCanvas: true, flipY: true });

            this.depthPickTarget.bind();
            this.renderVariant(renderer, camera, scene, helper, 'depth', PickType.None);
        }
    }
}

let AsyncPickingWarningShown = false;

export function checkAsyncPickingSupport(webgl: WebGLContext): boolean {
    if (webgl.isWebGL2) return true;

    if (isDebugMode && !AsyncPickingWarningShown) {
        console.log('WebGL2 required for async picking. Falling back to synchronous picking.');
        AsyncPickingWarningShown = true;
    }
    return false;
}

export enum AsyncPickStatus { Pending, Resolved, Failed };

export class PickBuffers {
    private object: Uint8Array;
    private instance: Uint8Array;
    private group: Uint8Array;
    private depth: Uint8Array;

    private objectBuffer: PixelPackBuffer;
    private instanceBuffer: PixelPackBuffer;
    private groupBuffer: PixelPackBuffer;
    private depthBuffer: PixelPackBuffer;

    private viewport = Viewport.create(0, 0, 0, 0);

    private setup() {
        const size = this.viewport.width * this.viewport.height * 4;
        if (!this.object || this.object.length !== size) {
            this.object = new Uint8Array(size);
            this.instance = new Uint8Array(size);
            this.group = new Uint8Array(size);
            this.depth = new Uint8Array(size);
        }
    }

    setViewport(x: number, y: number, width: number, height: number) {
        Viewport.set(this.viewport, x, y, width, height);
        this.setup();
    }

    read() {
        if (isTimingMode) this.webgl.timer.mark('PickBuffers.read');
        const { x, y, width, height } = this.viewport;

        this.pickPass.bindObject();
        this.webgl.readPixels(x, y, width, height, this.object);

        this.pickPass.bindInstance();
        this.webgl.readPixels(x, y, width, height, this.instance);

        this.pickPass.bindGroup();
        this.webgl.readPixels(x, y, width, height, this.group);

        this.pickPass.bindDepth();
        this.webgl.readPixels(x, y, width, height, this.depth);

        this.ready = true;
        if (isTimingMode) this.webgl.timer.markEnd('PickBuffers.read');
    }

    private fenceSync: WebGLSync | null = null;
    private fenceTimestamp: number = 0;

    private ready = false;
    private lag = 0;

    asyncRead() {
        const { gl } = this.webgl;
        if (!isWebGL2(gl)) return;

        if (isTimingMode) this.webgl.timer.mark('PickBuffers.asyncRead');
        if (this.fenceSync !== null) {
            gl.deleteSync(this.fenceSync);
        }
        const { x, y, width, height } = this.viewport;

        this.pickPass.bindObject();
        this.objectBuffer.read(x, y, width, height);

        this.pickPass.bindInstance();
        this.instanceBuffer.read(x, y, width, height);

        this.pickPass.bindGroup();
        this.groupBuffer.read(x, y, width, height);

        this.pickPass.bindDepth();
        this.depthBuffer.read(x, y, width, height);

        this.fenceTimestamp = now();
        this.fenceSync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
        // gl.flush();

        this.ready = false;
        if (isTimingMode) this.webgl.timer.markEnd('PickBuffers.asyncRead');
    }

    check(): AsyncPickStatus {
        if (this.ready) return AsyncPickStatus.Resolved;
        if (this.fenceSync === null) return AsyncPickStatus.Failed;

        const { gl } = this.webgl;
        if (!isWebGL2(gl)) return AsyncPickStatus.Failed;

        const res = gl.clientWaitSync(this.fenceSync, 0, 0);
        if (res === gl.WAIT_FAILED || this.lag >= this.maxAsyncReadLag) {
            // console.log(`failed to get buffer data after ${this.lag + 1} checks`);
            if (res !== gl.WAIT_FAILED && now() - this.fenceTimestamp < 1000 / 60) {
                this.lag += 1;
                return AsyncPickStatus.Pending;
            }
            gl.deleteSync(this.fenceSync);
            this.fenceSync = null;
            this.lag = 0;
            this.ready = false;
            return AsyncPickStatus.Failed;
        } else if (res === gl.TIMEOUT_EXPIRED) {
            this.lag += 1;
            // console.log(`waiting for buffer data for ${this.lag} checks`);
            return AsyncPickStatus.Pending;
        } else {
            this.objectBuffer.getSubData(this.object);
            this.instanceBuffer.getSubData(this.instance);
            this.groupBuffer.getSubData(this.group);
            this.depthBuffer.getSubData(this.depth);

            // console.log(`got buffer data after ${this.lag + 1} checks`);
            gl.deleteSync(this.fenceSync);
            this.fenceSync = null;
            this.lag = 0;
            this.ready = true;

            return AsyncPickStatus.Resolved;
        }
    }

    private getIdx(x: number, y: number): number {
        return (y * this.viewport.width + x) * 4;
    }

    getDepth(x: number, y: number): number {
        if (!this.ready) return -1;

        const idx = this.getIdx(x, y);
        const b = this.depth;
        return unpackRGBAToDepth(b[idx], b[idx + 1], b[idx + 2], b[idx + 3]);
    }

    private getId(x: number, y: number, buffer: Uint8Array) {
        if (!this.ready) return -1;

        const idx = this.getIdx(x, y);
        return unpackRGBToInt(buffer[idx], buffer[idx + 1], buffer[idx + 2]);
    }

    getObjectId(x: number, y: number) {
        return this.getId(x, y, this.object);
    }

    getInstanceId(x: number, y: number) {
        return this.getId(x, y, this.instance);
    }

    getGroupId(x: number, y: number) {
        return this.getId(x, y, this.group);
    }

    getPickingId(x: number, y: number): PickingId | undefined {
        const objectId = this.getObjectId(x, y);
        // console.log('objectId', objectId);
        if (objectId === -1 || objectId === NullId) return;

        const instanceId = this.getInstanceId(x, y);
        // console.log('instanceId', instanceId);
        if (instanceId === -1 || instanceId === NullId) return;

        const groupId = this.getGroupId(x, y);
        // console.log('groupId', groupId);
        if (groupId === -1 || groupId === NullId) return;

        return { objectId, instanceId, groupId };
    }

    reset() {
        this.fenceSync = null;
        this.ready = false;
        this.lag = 0;
        this.fenceTimestamp = 0;
    }

    dispose() {
        const { gl } = this.webgl;
        if (!isWebGL2(gl)) return;

        this.objectBuffer.destroy();
        this.instanceBuffer.destroy();
        this.groupBuffer.destroy();
        this.depthBuffer.destroy();

        if (this.fenceSync !== null) {
            gl.deleteSync(this.fenceSync);
            this.fenceSync = null;
        }
    }

    constructor(private webgl: WebGLContext, private pickPass: PickPass, public maxAsyncReadLag = 5) {
        if (webgl.isWebGL2) {
            this.objectBuffer = webgl.resources.pixelPack('rgba', 'ubyte');
            this.instanceBuffer = webgl.resources.pixelPack('rgba', 'ubyte');
            this.groupBuffer = webgl.resources.pixelPack('rgba', 'ubyte');
            this.depthBuffer = webgl.resources.pixelPack('rgba', 'ubyte');
        }

        this.setup();
    }
}
