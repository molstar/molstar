/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PickingId } from '../../mol-geo/geometry/picking';
import { Renderer } from '../../mol-gl/renderer';
import { Scene } from '../../mol-gl/scene';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { GraphicsRenderVariant } from '../../mol-gl/webgl/render-item';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { Vec3 } from '../../mol-math/linear-algebra';
import { spiral2d } from '../../mol-math/misc';
import { decodeFloatRGB, unpackRGBAToDepth } from '../../mol-util/float-packing';
import { Camera, ICamera } from '../camera';
import { StereoCamera } from '../camera/stereo';
import { cameraUnproject } from '../camera/util';
import { Viewport } from '../camera/util';
import { Helper } from '../helper/helper';
import { DrawPass } from './draw';

const NullId = Math.pow(2, 24) - 2;

export type PickData = { id: PickingId, position: Vec3 }

export class PickPass {
    readonly objectPickTarget: RenderTarget
    readonly instancePickTarget: RenderTarget
    readonly groupPickTarget: RenderTarget
    readonly depthPickTarget: RenderTarget

    private pickWidth: number
    private pickHeight: number

    constructor(private webgl: WebGLContext, private drawPass: DrawPass, readonly pickBaseScale: number) {
        const pickScale = pickBaseScale / webgl.pixelRatio;
        this.pickWidth = Math.ceil(drawPass.colorTarget.getWidth() * pickScale);
        this.pickHeight = Math.ceil(drawPass.colorTarget.getHeight() * pickScale);

        this.objectPickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
        this.instancePickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
        this.groupPickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
        this.depthPickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
    }

    get drawingBufferHeight() {
        return this.drawPass.colorTarget.getHeight();
    }

    syncSize() {
        const pickScale = this.pickBaseScale / this.webgl.pixelRatio;
        const pickWidth = Math.ceil(this.drawPass.colorTarget.getWidth() * pickScale);
        const pickHeight = Math.ceil(this.drawPass.colorTarget.getHeight() * pickScale);

        if (pickWidth !== this.pickWidth || pickHeight !== this.pickHeight) {
            this.pickWidth = pickWidth;
            this.pickHeight = pickHeight;

            this.objectPickTarget.setSize(this.pickWidth, this.pickHeight);
            this.instancePickTarget.setSize(this.pickWidth, this.pickHeight);
            this.groupPickTarget.setSize(this.pickWidth, this.pickHeight);
            this.depthPickTarget.setSize(this.pickWidth, this.pickHeight);
        }
    }

    private renderVariant(renderer: Renderer, camera: ICamera, scene: Scene, helper: Helper, variant: GraphicsRenderVariant) {
        const depth = this.drawPass.depthTexturePrimitives;
        renderer.clear(false);

        renderer.update(camera);
        renderer.renderPick(scene.primitives, camera, variant, null);
        renderer.renderPick(scene.volumes, camera, variant, depth);
        renderer.renderPick(helper.handle.scene, camera, variant, null);

        if (helper.camera.isEnabled) {
            helper.camera.update(camera);
            renderer.update(helper.camera.camera);
            renderer.renderPick(helper.camera.scene, helper.camera.camera, variant, null);
        }
    }

    render(renderer: Renderer, camera: ICamera, scene: Scene, helper: Helper) {
        this.objectPickTarget.bind();
        this.renderVariant(renderer, camera, scene, helper, 'pickObject');

        this.instancePickTarget.bind();
        this.renderVariant(renderer, camera, scene, helper, 'pickInstance');

        this.groupPickTarget.bind();
        this.renderVariant(renderer, camera, scene, helper, 'pickGroup');
        // printTexture(this.webgl, this.groupPickTarget.texture, { id: 'group' })

        this.depthPickTarget.bind();
        this.renderVariant(renderer, camera, scene, helper, 'depth');
    }
}

export class PickHelper {
    dirty = true

    private objectBuffer: Uint8Array
    private instanceBuffer: Uint8Array
    private groupBuffer: Uint8Array
    private depthBuffer: Uint8Array

    private viewport = Viewport()

    private pickScale: number
    private pickX: number
    private pickY: number
    private pickWidth: number
    private pickHeight: number
    private halfPickWidth: number

    private spiral: [number, number][]

    private setupBuffers() {
        const bufferSize = this.pickWidth * this.pickHeight * 4;
        if (!this.objectBuffer || this.objectBuffer.length !== bufferSize) {
            this.objectBuffer = new Uint8Array(bufferSize);
            this.instanceBuffer = new Uint8Array(bufferSize);
            this.groupBuffer = new Uint8Array(bufferSize);
            this.depthBuffer = new Uint8Array(bufferSize);
        }
    }

    setViewport(x: number, y: number, width: number, height: number) {
        Viewport.set(this.viewport, x, y, width, height);

        this.pickScale = this.pickPass.pickBaseScale / this.webgl.pixelRatio;
        this.pickX = Math.ceil(x * this.pickScale);
        this.pickY = Math.ceil(y * this.pickScale);

        const pickWidth = Math.floor(width * this.pickScale);
        const pickHeight = Math.floor(height * this.pickScale);

        if (pickWidth !== this.pickWidth || pickHeight !== this.pickHeight) {
            this.pickWidth = pickWidth;
            this.pickHeight = pickHeight;
            this.halfPickWidth = Math.floor(this.pickWidth / 2);

            this.setupBuffers();
        }

        this.spiral = spiral2d(Math.round(this.pickScale * this.pickPadding));
    }

    private syncBuffers() {
        const { pickX, pickY, pickWidth, pickHeight } = this;

        this.pickPass.objectPickTarget.bind();
        this.webgl.readPixels(pickX, pickY, pickWidth, pickHeight, this.objectBuffer);

        this.pickPass.instancePickTarget.bind();
        this.webgl.readPixels(pickX, pickY, pickWidth, pickHeight, this.instanceBuffer);

        this.pickPass.groupPickTarget.bind();
        this.webgl.readPixels(pickX, pickY, pickWidth, pickHeight, this.groupBuffer);

        this.pickPass.depthPickTarget.bind();
        this.webgl.readPixels(pickX, pickY, pickWidth, pickHeight, this.depthBuffer);
    }

    private getBufferIdx(x: number, y: number): number {
        return (y * this.pickWidth + x) * 4;
    }

    private getDepth(x: number, y: number): number {
        const idx = this.getBufferIdx(x, y);
        const b = this.depthBuffer;
        return unpackRGBAToDepth(b[idx], b[idx + 1], b[idx + 2], b[idx + 3]);
    }

    private getId(x: number, y: number, buffer: Uint8Array) {
        const idx = this.getBufferIdx(x, y);
        return decodeFloatRGB(buffer[idx], buffer[idx + 1], buffer[idx + 2]);
    }

    private render(camera: Camera | StereoCamera) {
        const { pickX, pickY, pickWidth, pickHeight, halfPickWidth } = this;
        const { renderer, scene, helper } = this;

        renderer.setTransparentBackground(false);
        renderer.setDrawingBufferSize(this.pickPass.objectPickTarget.getWidth(), this.pickPass.objectPickTarget.getHeight());
        renderer.setPixelRatio(this.pickScale);

        if (StereoCamera.is(camera)) {
            renderer.setViewport(pickX, pickY, halfPickWidth, pickHeight);
            this.pickPass.render(renderer, camera.left, scene, helper);

            renderer.setViewport(pickX + halfPickWidth, pickY, pickWidth - halfPickWidth, pickHeight);
            this.pickPass.render(renderer, camera.right, scene, helper);
        } else {
            renderer.setViewport(pickX, pickY, pickWidth, pickHeight);
            this.pickPass.render(renderer, camera, scene, helper);
        }

        this.dirty = false;
    }

    private identifyInternal(x: number, y: number, camera: Camera | StereoCamera): PickData | undefined {
        const { webgl, pickScale } = this;
        if (webgl.isContextLost) return;

        x *= webgl.pixelRatio;
        y *= webgl.pixelRatio;
        y = this.pickPass.drawingBufferHeight - y; // flip y

        const { viewport } = this;

        // check if within viewport
        if (x < viewport.x ||
            y < viewport.y ||
            x > viewport.x + viewport.width ||
            y > viewport.y + viewport.height
        ) return;

        if (this.dirty) {
            this.render(camera);
            this.syncBuffers();
        }

        const xv = x - viewport.x;
        const yv = y - viewport.y;

        const xp = Math.floor(xv * pickScale);
        const yp = Math.floor(yv * pickScale);

        const objectId = this.getId(xp, yp, this.objectBuffer);
        // console.log('objectId', objectId);
        if (objectId === -1 || objectId === NullId) return;

        const instanceId = this.getId(xp, yp, this.instanceBuffer);
        // console.log('instanceId', instanceId);
        if (instanceId === -1 || instanceId === NullId) return;

        const groupId = this.getId(xp, yp, this.groupBuffer);
        // console.log('groupId', groupId);
        if (groupId === -1 || groupId === NullId) return;

        const z = this.getDepth(xp, yp);
        const position = Vec3.create(x, viewport.height - y, z);
        if (StereoCamera.is(camera)) {
            const halfWidth = Math.floor(viewport.width / 2);
            if (x > viewport.x + halfWidth) {
                position[0] = viewport.x + (xv - halfWidth) * 2;
                cameraUnproject(position, position, viewport, camera.right.inverseProjectionView);
            } else {
                position[0] = viewport.x + xv * 2;
                cameraUnproject(position, position, viewport, camera.left.inverseProjectionView);
            }
        } else {
            cameraUnproject(position, position, viewport, camera.inverseProjectionView);
        }

        // console.log({ { objectId, instanceId, groupId }, position} );
        return { id: { objectId, instanceId, groupId }, position };
    }

    identify(x: number, y: number, camera: Camera | StereoCamera): PickData | undefined {
        for (const d of this.spiral) {
            const pickData = this.identifyInternal(x + d[0], y + d[1], camera);
            if (pickData) return pickData;
        }
    }

    constructor(private webgl: WebGLContext, private renderer: Renderer, private scene: Scene, private helper: Helper, private pickPass: PickPass, viewport: Viewport, readonly pickPadding = 1) {
        this.setViewport(viewport.x, viewport.y, viewport.width, viewport.height);
    }
}