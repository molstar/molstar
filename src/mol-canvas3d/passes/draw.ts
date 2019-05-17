/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from 'mol-gl/webgl/context';
import { createRenderTarget, RenderTarget } from 'mol-gl/webgl/render-target';
import Renderer from 'mol-gl/renderer';
import Scene from 'mol-gl/scene';
import { BoundingSphereHelper } from '../helper/bounding-sphere-helper';
import { createTexture, Texture } from 'mol-gl/webgl/texture';


export class DrawPass {
    colorTarget: RenderTarget
    depthTexture: Texture
    packedDepth: boolean

    private depthTarget: RenderTarget | null

    constructor(private webgl: WebGLContext, private renderer: Renderer, private scene: Scene, private debugHelper: BoundingSphereHelper) {
        const { gl, extensions } = webgl
        const width = gl.drawingBufferWidth
        const height = gl.drawingBufferHeight
        this.colorTarget = createRenderTarget(webgl, gl.drawingBufferWidth, gl.drawingBufferHeight)
        this.packedDepth = !extensions.depthTexture
        this.depthTarget = this.packedDepth ? createRenderTarget(webgl, width, height) : null
        this.depthTexture = this.depthTarget ? this.depthTarget.texture : createTexture(webgl, 'image-depth', 'depth', 'ushort', 'nearest')
        if (!this.packedDepth) {
            this.depthTexture.define(width, height)
            this.depthTexture.attachFramebuffer(this.colorTarget.framebuffer, 'depth')
        }
    }

    setSize(width: number, height: number) {
        this.colorTarget.setSize(width, height)
        if (this.depthTarget) {
            this.depthTarget.setSize(width, height)
        } else {
            this.depthTexture.define(width, height)
        }
    }

    render(toDrawingBuffer: boolean) {
        const { webgl, renderer, scene, debugHelper, colorTarget, depthTarget } = this
        const { gl } = webgl
        if (toDrawingBuffer) {
            webgl.unbindFramebuffer()
            gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight)
        } else {
            colorTarget.bind()
        }
        renderer.setViewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight)
        renderer.render(scene, 'color', true)
        if (debugHelper.isEnabled) {
            debugHelper.syncVisibility()
            renderer.render(debugHelper.scene, 'color', false)
        }

        // do a depth pass if not rendering to drawing buffer and
        // extensions.depthTexture is unsupported (i.e. depthTarget is set)
        if (!toDrawingBuffer && depthTarget) {
            depthTarget.bind()
            renderer.render(scene, 'depth', true)
            if (debugHelper.isEnabled) {
                debugHelper.syncVisibility()
                renderer.render(debugHelper.scene, 'depth', false)
            }
        }
    }
}