/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PickingId } from '../../mol-geo/geometry/picking';
import Renderer from '../../mol-gl/renderer';
import Scene from '../../mol-gl/scene';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { GraphicsRenderVariant } from '../../mol-gl/webgl/render-item';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { Vec3 } from '../../mol-math/linear-algebra';
import { decodeFloatRGB, unpackRGBAToDepth } from '../../mol-util/float-packing';
import { Camera, ICamera } from '../camera';
import { StereoCamera } from '../camera/stereo';
import { cameraUnproject } from '../camera/util';
import { HandleHelper } from '../helper/handle-helper';
import { DrawPass } from './draw';

const NullId = Math.pow(2, 24) - 2;

export type PickData = { id: PickingId, position: Vec3 }

export class PickPass {
    pickDirty = true

    objectPickTarget: RenderTarget
    instancePickTarget: RenderTarget
    groupPickTarget: RenderTarget
    depthPickTarget: RenderTarget

    isStereo = false

    private objectBuffer: Uint8Array
    private instanceBuffer: Uint8Array
    private groupBuffer: Uint8Array
    private depthBuffer: Uint8Array

    private pickScale: number
    private pickWidth: number
    private pickHeight: number

    constructor(private webgl: WebGLContext, private renderer: Renderer, private scene: Scene, private camera: Camera, private stereoCamera: StereoCamera, private handleHelper: HandleHelper, private pickBaseScale: number, private drawPass: DrawPass) {
        this.pickScale = pickBaseScale / webgl.pixelRatio;
        this.pickWidth = Math.ceil(camera.viewport.width * this.pickScale);
        this.pickHeight = Math.ceil(camera.viewport.height * this.pickScale);

        this.objectPickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
        this.instancePickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
        this.groupPickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
        this.depthPickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);

        this.setupBuffers();
    }

    private setupBuffers() {
        const bufferSize = this.pickWidth * this.pickHeight * 4;
        if (!this.objectBuffer || this.objectBuffer.length !== bufferSize) {
            this.objectBuffer = new Uint8Array(bufferSize);
            this.instanceBuffer = new Uint8Array(bufferSize);
            this.groupBuffer = new Uint8Array(bufferSize);
            this.depthBuffer = new Uint8Array(bufferSize);
        }
    }

    setSize(width: number, height: number) {
        this.pickScale = this.pickBaseScale / this.webgl.pixelRatio;
        const pickWidth = Math.ceil(width * this.pickScale);
        const pickHeight = Math.ceil(height * this.pickScale);

        if (pickWidth !== this.pickWidth || pickHeight !== this.pickHeight) {
            this.pickWidth = pickWidth;
            this.pickHeight = pickHeight;

            this.objectPickTarget.setSize(this.pickWidth, this.pickHeight);
            this.instancePickTarget.setSize(this.pickWidth, this.pickHeight);
            this.groupPickTarget.setSize(this.pickWidth, this.pickHeight);
            this.depthPickTarget.setSize(this.pickWidth, this.pickHeight);

            this.setupBuffers();
        }
    }

    private renderVariant(variant: GraphicsRenderVariant) {
        if (this.isStereo) {
            const w = (this.pickWidth / 2) | 0;

            this.renderer.setViewport(0, 0, w, this.pickHeight);
            this._renderVariant(this.stereoCamera.left, variant);

            this.renderer.setViewport(w, 0, this.pickWidth - w, this.pickHeight);
            this._renderVariant(this.stereoCamera.right, variant);
        } else {
            this.renderer.setViewport(0, 0, this.pickWidth, this.pickHeight);
            this._renderVariant(this.camera, variant);
        }
    }

    private _renderVariant(camera: ICamera, variant: GraphicsRenderVariant) {
        const { renderer, scene, handleHelper: { scene: handleScene } } = this;
        const depth = this.drawPass.depthTexturePrimitives;

        renderer.render(scene.primitives, camera, variant, true, false, null);
        renderer.render(scene.volumes, camera, variant, false, false, depth);
        renderer.render(handleScene, camera, variant, false, false, null);
    }

    render() {
        this.objectPickTarget.bind();
        this.renderVariant('pickObject');

        this.instancePickTarget.bind();
        this.renderVariant('pickInstance');

        this.groupPickTarget.bind();
        this.renderVariant('pickGroup');

        this.depthPickTarget.bind();
        this.renderVariant('depth');

        this.pickDirty = false;
    }

    private syncBuffers() {
        const { webgl } = this;

        this.objectPickTarget.bind();
        webgl.readPixels(0, 0, this.pickWidth, this.pickHeight, this.objectBuffer);

        this.instancePickTarget.bind();
        webgl.readPixels(0, 0, this.pickWidth, this.pickHeight, this.instanceBuffer);

        this.groupPickTarget.bind();
        webgl.readPixels(0, 0, this.pickWidth, this.pickHeight, this.groupBuffer);

        this.depthPickTarget.bind();
        webgl.readPixels(0, 0, this.pickWidth, this.pickHeight, this.depthBuffer);
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

    identify(x: number, y: number): PickData | undefined {
        const { webgl, pickScale, camera: { viewport } } = this;
        if (webgl.isContextLost) return;

        const { gl, pixelRatio } = webgl;
        x *= pixelRatio;
        y *= pixelRatio;

        // check if within viewport
        if (x < viewport.x ||
            gl.drawingBufferHeight - y < viewport.y ||
            x > viewport.x + viewport.width ||
            gl.drawingBufferHeight - y > viewport.y + viewport.height
        ) {
            return;
        }

        if (this.pickDirty) {
            this.render();
            this.syncBuffers();
        }

        x -= viewport.x;
        y += viewport.y; // plus because of flipped y
        y = gl.drawingBufferHeight - y; // flip y

        const xp = Math.floor(x * pickScale);
        const yp = Math.floor(y * pickScale);

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
        const position = Vec3.create(x, gl.drawingBufferHeight - y, z);
        cameraUnproject(position, position, viewport, this.camera.inverseProjectionView);

        // console.log({ { objectId, instanceId, groupId }, position} );
        return { id: { objectId, instanceId, groupId }, position };
    }
}