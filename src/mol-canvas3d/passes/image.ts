/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { Renderer } from '../../mol-gl/renderer';
import { Scene } from '../../mol-gl/scene';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { DrawPass } from './draw';
import { PostprocessingParams } from './postprocessing';
import { MultiSamplePass, MultiSampleParams, MultiSampleHelper } from './multi-sample';
import { Camera } from '../camera';
import { Viewport } from '../camera/util';
import { PixelData } from '../../mol-util/image';
import { Helper } from '../helper/helper';
import { CameraHelper, CameraHelperParams } from '../helper/camera-helper';
import { MarkingParams } from './marking';

export const ImageParams = {
    transparentBackground: PD.Boolean(false),
    multiSample: PD.Group(MultiSampleParams),
    postprocessing: PD.Group(PostprocessingParams),
    marking: PD.Group(MarkingParams),

    cameraHelper: PD.Group(CameraHelperParams),
};
export type ImageProps = PD.Values<typeof ImageParams>

export class ImagePass {
    private _width = 0
    private _height = 0
    private _camera = new Camera()

    readonly props: ImageProps

    private _colorTarget: RenderTarget
    get colorTarget() { return this._colorTarget; }

    private readonly drawPass: DrawPass
    private readonly multiSamplePass: MultiSamplePass
    private readonly multiSampleHelper: MultiSampleHelper
    private readonly helper: Helper

    get width() { return this._width; }
    get height() { return this._height; }

    constructor(private webgl: WebGLContext, private renderer: Renderer, private scene: Scene, private camera: Camera, helper: Helper, enableWboit: boolean, props: Partial<ImageProps>) {
        this.props = { ...PD.getDefaultValues(ImageParams), ...props };

        this.drawPass = new DrawPass(webgl, 128, 128, enableWboit);
        this.multiSamplePass = new MultiSamplePass(webgl, this.drawPass);
        this.multiSampleHelper = new MultiSampleHelper(this.multiSamplePass);

        this.helper = {
            camera: new CameraHelper(webgl, this.props.cameraHelper),
            debug: helper.debug,
            handle: helper.handle,
        };

        this.setSize(1024, 768);
    }

    setSize(width: number, height: number) {
        if (width === this._width && height === this._height) return;

        this._width = width;
        this._height = height;

        this.drawPass.setSize(width, height);
        this.multiSamplePass.syncSize();
    }

    setProps(props: Partial<ImageProps> = {}) {
        Object.assign(this.props, props);
        if (props.cameraHelper) this.helper.camera.setProps(props.cameraHelper);
    }

    render() {
        Camera.copySnapshot(this._camera.state, this.camera.state);
        Viewport.set(this._camera.viewport, 0, 0, this._width, this._height);
        this._camera.update();

        if (MultiSamplePass.isEnabled(this.props.multiSample)) {
            this.multiSampleHelper.render(this.renderer, this._camera, this.scene, this.helper, false, this.props.transparentBackground, this.props);
            this._colorTarget = this.multiSamplePass.colorTarget;
        } else {
            this.drawPass.render(this.renderer, this._camera, this.scene, this.helper, false, this.props.transparentBackground, this.props.postprocessing, this.props.marking);
            this._colorTarget = this.drawPass.getColorTarget(this.props.postprocessing);
        }
    }

    getImageData(width: number, height: number, viewport?: Viewport) {
        this.setSize(width, height);
        this.render();
        this.colorTarget.bind();

        const w = viewport?.width ?? width, h = viewport?.height ?? height;

        const array = new Uint8Array(w * h * 4);
        if (!viewport) {
            this.webgl.readPixels(0, 0, w, h, array);
        } else {
            this.webgl.readPixels(viewport.x, height - viewport.y - viewport.height, w, h, array);
        }
        const pixelData = PixelData.create(array, w, h);
        PixelData.flipY(pixelData);
        PixelData.divideByAlpha(pixelData);
        return new ImageData(new Uint8ClampedArray(array), w, h);
    }
}