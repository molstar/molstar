/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from 'mol-gl/webgl/context';
import { createRenderTarget, RenderTarget } from 'mol-gl/webgl/render-target';
import Renderer from 'mol-gl/renderer';
import Scene from 'mol-gl/scene';
import { PickingId } from 'mol-geo/geometry/picking';
import { decodeFloatRGB } from 'mol-util/float-packing';

const readBuffer = new Uint8Array(4)

export class PickPass {
    pickDirty = true
    objectPickTarget: RenderTarget
    instancePickTarget: RenderTarget
    groupPickTarget: RenderTarget

    private pickScale: number
    private pickWidth: number
    private pickHeight: number

    constructor(private webgl: WebGLContext, private renderer: Renderer, private scene: Scene, private pickBaseScale: number) {
        const { gl } = webgl
        const width = gl.drawingBufferWidth
        const height = gl.drawingBufferHeight

        this.pickScale = pickBaseScale / webgl.pixelRatio
        this.pickWidth = Math.round(width * this.pickScale)
        this.pickHeight = Math.round(height * this.pickScale)
        this.objectPickTarget = createRenderTarget(webgl, this.pickWidth, this.pickHeight)
        this.instancePickTarget = createRenderTarget(webgl, this.pickWidth, this.pickHeight)
        this.groupPickTarget = createRenderTarget(webgl, this.pickWidth, this.pickHeight)
    }

    setSize(width: number, height: number) {
        this.pickScale = this.pickBaseScale / this.webgl.pixelRatio
        this.pickWidth = Math.round(width * this.pickScale)
        this.pickHeight = Math.round(height * this.pickScale)
        this.objectPickTarget.setSize(this.pickWidth, this.pickHeight)
        this.instancePickTarget.setSize(this.pickWidth, this.pickHeight)
        this.groupPickTarget.setSize(this.pickWidth, this.pickHeight)
    }

    render() {
        const { renderer, scene } = this
        renderer.setViewport(0, 0, this.pickWidth, this.pickHeight);
        this.objectPickTarget.bind();
        renderer.render(scene, 'pickObject', true);
        this.instancePickTarget.bind();
        renderer.render(scene, 'pickInstance', true);
        this.groupPickTarget.bind();
        renderer.render(scene, 'pickGroup', true);
    }

    identify(x: number, y: number): PickingId | undefined {
        const { webgl, pickScale } = this
        const { gl } = webgl
        if (this.pickDirty) this.render()

        x *= webgl.pixelRatio
        y *= webgl.pixelRatio
        y = gl.drawingBufferHeight - y // flip y

        const xp = Math.round(x * pickScale)
        const yp = Math.round(y * pickScale)

        this.objectPickTarget.bind()
        webgl.readPixels(xp, yp, 1, 1, readBuffer)
        const objectId = decodeFloatRGB(readBuffer[0], readBuffer[1], readBuffer[2])
        if (objectId === -1) return

        this.instancePickTarget.bind()
        webgl.readPixels(xp, yp, 1, 1, readBuffer)
        const instanceId = decodeFloatRGB(readBuffer[0], readBuffer[1], readBuffer[2])
        if (instanceId === -1) return

        this.groupPickTarget.bind()
        webgl.readPixels(xp, yp, 1, 1, readBuffer)
        const groupId = decodeFloatRGB(readBuffer[0], readBuffer[1], readBuffer[2])
        if (groupId === -1) return

        return { objectId, instanceId, groupId }
    }
}