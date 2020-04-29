/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import Renderer from '../../mol-gl/renderer';
import Scene from '../../mol-gl/scene';
import { BoundingSphereHelper } from '../helper/bounding-sphere-helper';
import { Texture } from '../../mol-gl/webgl/texture';
import { Camera } from '../camera';
import { CameraHelper, CameraHelperParams } from '../helper/camera-helper';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { HandleHelper } from '../helper/handle-helper';

export const DrawPassParams = {
    cameraHelper: PD.Group(CameraHelperParams)
};
export const DefaultDrawPassProps = PD.getDefaultValues(DrawPassParams);
export type DrawPassProps = PD.Values<typeof DrawPassParams>

export class DrawPass {
    colorTarget: RenderTarget
    depthTexture: Texture
    packedDepth: boolean

    cameraHelper: CameraHelper

    private depthTarget: RenderTarget | null

    constructor(private webgl: WebGLContext, private renderer: Renderer, private scene: Scene, private camera: Camera, private debugHelper: BoundingSphereHelper, private handleHelper: HandleHelper, props: Partial<DrawPassProps> = {}) {
        const { gl, extensions, resources } = webgl;
        const width = gl.drawingBufferWidth;
        const height = gl.drawingBufferHeight;
        this.colorTarget = webgl.createRenderTarget(width, height);
        this.packedDepth = !extensions.depthTexture;
        this.depthTarget = this.packedDepth ? webgl.createRenderTarget(width, height) : null;
        this.depthTexture = this.depthTarget ? this.depthTarget.texture : resources.texture('image-depth', 'depth', 'ushort', 'nearest');
        if (!this.packedDepth) {
            this.depthTexture.define(width, height);
            this.depthTexture.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
        }

        const p = { ...DefaultDrawPassProps, ...props };
        this.cameraHelper = new CameraHelper(webgl, p.cameraHelper);
    }

    setSize(width: number, height: number) {
        this.colorTarget.setSize(width, height);
        if (this.depthTarget) {
            this.depthTarget.setSize(width, height);
        } else {
            this.depthTexture.define(width, height);
        }
    }

    setProps(props: Partial<DrawPassProps>) {
        if (props.cameraHelper) this.cameraHelper.setProps(props.cameraHelper);
    }

    get props(): DrawPassProps {
        return {
            cameraHelper: { ...this.cameraHelper.props }
        };
    }

    render(toDrawingBuffer: boolean, transparentBackground: boolean) {
        const { webgl, renderer, colorTarget, depthTarget } = this;
        if (toDrawingBuffer) {
            webgl.unbindFramebuffer();
        } else {
            colorTarget.bind();
            if (!this.packedDepth) {
                // TODO unlcear why it is not enough to call `attachFramebuffer` in `Texture.reset`
                this.depthTexture.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
            }
        }

        renderer.setViewport(0, 0, colorTarget.getWidth(), colorTarget.getHeight());
        this.renderInternal('color', transparentBackground);

        // do a depth pass if not rendering to drawing buffer and
        // extensions.depthTexture is unsupported (i.e. depthTarget is set)
        if (!toDrawingBuffer && depthTarget) {
            depthTarget.bind();
            this.renderInternal('depth', transparentBackground);
        }
    }

    private renderInternal(variant: 'color' | 'depth', transparentBackground: boolean) {
        const { renderer, scene, camera, debugHelper, cameraHelper, handleHelper } = this;
        renderer.render(scene, camera, variant, true, transparentBackground);
        if (debugHelper.isEnabled) {
            debugHelper.syncVisibility();
            renderer.render(debugHelper.scene, camera, variant, false, transparentBackground);
        }
        if (handleHelper.isEnabled) {
            renderer.render(handleHelper.scene, camera, variant, false, transparentBackground);
        }
        if (cameraHelper.isEnabled) {
            cameraHelper.update(camera);
            renderer.render(cameraHelper.scene, cameraHelper.camera, variant, false, transparentBackground);
        }
    }
}