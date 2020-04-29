/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import Renderer from '../../mol-gl/renderer';
import Scene from '../../mol-gl/scene';
import { PickingId } from '../../mol-geo/geometry/picking';
import { decodeFloatRGB } from '../../mol-util/float-packing';
import { Camera } from '../camera';
import { HandleHelper } from '../helper/handle-helper';

export class PickPass {
    pickDirty = true

    objectPickTarget: RenderTarget
    instancePickTarget: RenderTarget
    groupPickTarget: RenderTarget

    private objectBuffer: Uint8Array
    private instanceBuffer: Uint8Array
    private groupBuffer: Uint8Array

    private pickScale: number
    private pickWidth: number
    private pickHeight: number

    constructor(private webgl: WebGLContext, private renderer: Renderer, private scene: Scene, private camera: Camera, private handleHelper: HandleHelper, private pickBaseScale: number) {
        const { gl } = webgl;
        const width = gl.drawingBufferWidth;
        const height = gl.drawingBufferHeight;

        this.pickScale = pickBaseScale / webgl.pixelRatio;
        this.pickWidth = Math.round(width * this.pickScale);
        this.pickHeight = Math.round(height * this.pickScale);

        this.objectPickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
        this.instancePickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);
        this.groupPickTarget = webgl.createRenderTarget(this.pickWidth, this.pickHeight);

        this.setupBuffers();
    }

    private setupBuffers() {
        const bufferSize = this.pickWidth * this.pickHeight * 4;
        if (!this.objectBuffer || this.objectBuffer.length !== bufferSize) {
            this.objectBuffer = new Uint8Array(bufferSize);
            this.instanceBuffer = new Uint8Array(bufferSize);
            this.groupBuffer = new Uint8Array(bufferSize);
        }
    }

    setSize(width: number, height: number) {
        this.pickScale = this.pickBaseScale / this.webgl.pixelRatio;
        this.pickWidth = Math.round(width * this.pickScale);
        this.pickHeight = Math.round(height * this.pickScale);

        this.objectPickTarget.setSize(this.pickWidth, this.pickHeight);
        this.instancePickTarget.setSize(this.pickWidth, this.pickHeight);
        this.groupPickTarget.setSize(this.pickWidth, this.pickHeight);

        this.setupBuffers();
    }

    render() {
        const { renderer, scene, camera, handleHelper: { scene: handleScene } } = this;
        renderer.setViewport(0, 0, this.pickWidth, this.pickHeight);

        this.objectPickTarget.bind();
        renderer.render(scene, camera, 'pickObject', true, false);
        renderer.render(handleScene, camera, 'pickObject', false, false);
        this.instancePickTarget.bind();
        renderer.render(scene, camera, 'pickInstance', true, false);
        renderer.render(handleScene, camera, 'pickInstance', false, false);
        this.groupPickTarget.bind();
        renderer.render(scene, camera, 'pickGroup', true, false);
        renderer.render(handleScene, camera, 'pickGroup', false, false);

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
    }

    private getId(x: number, y: number, buffer: Uint8Array) {
        const idx = (y * this.pickWidth + x) * 4;
        return decodeFloatRGB(buffer[idx], buffer[idx + 1], buffer[idx + 2]);
    }

    identify(x: number, y: number): PickingId | undefined {
        const { webgl, pickScale } = this;
        if (webgl.isContextLost) return;

        const { gl } = webgl;
        if (this.pickDirty) {
            this.render();
            this.syncBuffers();
        }

        x *= webgl.pixelRatio;
        y *= webgl.pixelRatio;
        y = gl.drawingBufferHeight - y; // flip y

        const xp = Math.round(x * pickScale);
        const yp = Math.round(y * pickScale);

        const objectId = this.getId(xp, yp, this.objectBuffer);
        if (objectId === -1) return;

        const instanceId = this.getId(xp, yp, this.instanceBuffer);
        if (instanceId === -1) return;

        const groupId = this.getId(xp, yp, this.groupBuffer);
        if (groupId === -1) return;

        return { objectId, instanceId, groupId };
    }
}