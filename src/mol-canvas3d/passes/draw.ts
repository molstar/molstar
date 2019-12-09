/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { createRenderTarget, RenderTarget } from '../../mol-gl/webgl/render-target';
import Renderer from '../../mol-gl/renderer';
import Scene from '../../mol-gl/scene';
import { BoundingSphereHelper } from '../helper/bounding-sphere-helper';
import { createTexture, Texture } from '../../mol-gl/webgl/texture';
import { Camera } from '../camera';

export class DrawPass {
    colorTarget: RenderTarget
    depthTexture: Texture
    packedDepth: boolean

    private depthTarget: RenderTarget | null

    constructor(private webgl: WebGLContext, private renderer: Renderer, private scene: Scene, private camera: Camera, private debugHelper: BoundingSphereHelper) {
        const { gl, extensions } = webgl
        const width = gl.drawingBufferWidth
        const height = gl.drawingBufferHeight
        this.colorTarget = createRenderTarget(webgl, width, height)
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

    render(toDrawingBuffer: boolean, transparentBackground: boolean) {
        const { webgl, renderer, scene, camera, debugHelper, colorTarget, depthTarget } = this
        if (toDrawingBuffer) {
            webgl.unbindFramebuffer()
        } else {
            colorTarget.bind()
        }

        renderer.setViewport(0, 0, colorTarget.width, colorTarget.height)
        renderer.render(scene, camera, 'color', true, transparentBackground)
        if (debugHelper.isEnabled) {
            debugHelper.syncVisibility()
            renderer.render(debugHelper.scene, camera, 'color', false, transparentBackground)
        }

        // do a depth pass if not rendering to drawing buffer and
        // extensions.depthTexture is unsupported (i.e. depthTarget is set)
        if (!toDrawingBuffer && depthTarget) {
            depthTarget.bind()
            renderer.render(scene, camera, 'depth', true, transparentBackground)
            if (debugHelper.isEnabled) {
                debugHelper.syncVisibility()
                renderer.render(debugHelper.scene, camera, 'depth', false, transparentBackground)
            }
        }
    }
}