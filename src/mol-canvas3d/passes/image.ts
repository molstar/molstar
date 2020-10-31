/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import Renderer from '../../mol-gl/renderer';
import Scene from '../../mol-gl/scene';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { DrawPass } from './draw';
import { PostprocessingPass, PostprocessingParams } from './postprocessing';
import { MultiSamplePass, MultiSampleParams } from './multi-sample';
import { Camera } from '../camera';
import { Viewport } from '../camera/util';
import { PixelData } from '../../mol-util/image';
import { Helper } from '../helper/helper';
import { CameraHelper, CameraHelperParams } from '../helper/camera-helper';

export const ImageParams = {
    transparentBackground: PD.Boolean(false),
    multiSample: PD.Group(MultiSampleParams),
    postprocessing: PD.Group(PostprocessingParams),

    cameraHelper: PD.Group(CameraHelperParams),
};
export type ImageProps = PD.Values<typeof ImageParams>

export class ImagePass {
    private _width = 1024
    private _height = 768
    private _camera = new Camera()

    readonly props: ImageProps

    private _colorTarget: RenderTarget
    get colorTarget() { return this._colorTarget; }

    readonly drawPass: DrawPass

    private readonly postprocessing: PostprocessingPass
    private readonly multiSample: MultiSamplePass
    private readonly helper: Helper

    get width() { return this._width; }
    get height() { return this._height; }

    constructor(private webgl: WebGLContext, private renderer: Renderer, private scene: Scene, private camera: Camera, helper: Helper, props: Partial<ImageProps>) {
        this.props = { ...PD.getDefaultValues(ImageParams), ...props };

        this.drawPass = new DrawPass(webgl, 128, 128);
        this.postprocessing = new PostprocessingPass(webgl, this.drawPass);
        this.multiSample = new MultiSamplePass(webgl, this.drawPass, this.postprocessing);

        this.helper = {
            camera: new CameraHelper(webgl, this.props.cameraHelper),
            debug: helper.debug,
            handle: helper.handle,
        };

        this.setSize(this._width, this._height);
    }

    setSize(width: number, height: number) {
        if (width === this._width && height === this._height) return;

        this._width = width;
        this._height = height;

        this.drawPass.setSize(width, height);
        this.postprocessing.syncSize();
        this.multiSample.syncSize();
    }

    setProps(props: Partial<ImageProps> = {}) {
        Object.assign(this.props, props);
        if (props.cameraHelper) this.helper.camera.setProps(props.cameraHelper);
    }

    render() {
        Camera.copySnapshot(this._camera.state, this.camera.state);
        Viewport.set(this._camera.viewport, 0, 0, this._width, this._height);
        this._camera.update();

        this.renderer.setViewport(0, 0, this._width, this._height);

        if (MultiSamplePass.isEnabled(this.props.multiSample)) {
            this.multiSample.render(this.renderer, this._camera, this.scene, this.helper, false, this.props.transparentBackground, this.props);
            this._colorTarget = this.multiSample.colorTarget;
        } else {
            this.drawPass.render(this.renderer, this._camera, this.scene, this.helper, false, this.props.transparentBackground);
            if (PostprocessingPass.isEnabled(this.props.postprocessing)) {
                this.postprocessing.render(this._camera, false, this.props.postprocessing);
                this._colorTarget = this.postprocessing.target;
            } else {
                this._colorTarget = this.drawPass.colorTarget;
            }
        }
    }

    getImageData(width: number, height: number) {
        this.setSize(width, height);
        this.render();
        this.colorTarget.bind();
        const array = new Uint8Array(width * height * 4);
        this.webgl.readPixels(0, 0, width, height, array);
        PixelData.flipY({ array, width, height });
        return new ImageData(new Uint8ClampedArray(array), width, height);
    }
}