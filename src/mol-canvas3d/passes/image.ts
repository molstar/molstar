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
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { DrawPass, DrawPassParams } from './draw';
import { PostprocessingPass, PostprocessingParams } from './postprocessing';
import { MultiSamplePass, MultiSampleParams } from './multi-sample';
import { Camera } from '../camera';
import { Viewport } from '../camera/util';
import { HandleHelper } from '../helper/handle-helper';

export const ImageParams = {
    transparentBackground: PD.Boolean(false),
    multiSample: PD.Group(MultiSampleParams),
    postprocessing: PD.Group(PostprocessingParams),
    drawPass: PD.Group(DrawPassParams),
};
export type ImageProps = PD.Values<typeof ImageParams>

export class ImagePass {
    private _width = 1024
    private _height = 768
    private _camera = new Camera()
    private _transparentBackground = false

    private _colorTarget: RenderTarget
    get colorTarget() { return this._colorTarget; }

    readonly drawPass: DrawPass
    private readonly postprocessing: PostprocessingPass
    private readonly multiSample: MultiSamplePass

    get width() { return this._width; }
    get height() { return this._height; }

    constructor(webgl: WebGLContext, private renderer: Renderer, scene: Scene, private camera: Camera, debugHelper: BoundingSphereHelper, handleHelper: HandleHelper, props: Partial<ImageProps>) {
        const p = { ...PD.getDefaultValues(ImageParams), ...props };

        this._transparentBackground = p.transparentBackground;

        this.drawPass = new DrawPass(webgl, renderer, scene, this._camera, debugHelper, handleHelper, p.drawPass);
        this.postprocessing = new PostprocessingPass(webgl, this._camera, this.drawPass, p.postprocessing);
        this.multiSample = new MultiSamplePass(webgl, this._camera, this.drawPass, this.postprocessing, p.multiSample);

        this.setSize(this._width, this._height);
    }

    setSize(width: number, height: number) {
        if (width === this._width && height === this._height) return;

        this._width = width;
        this._height = height;

        this.drawPass.setSize(width, height);
        this.postprocessing.setSize(width, height);
        this.multiSample.setSize(width, height);
    }

    setProps(props: Partial<ImageProps> = {}) {
        if (props.transparentBackground !== undefined) this._transparentBackground = props.transparentBackground;
        if (props.postprocessing) this.postprocessing.setProps(props.postprocessing);
        if (props.multiSample) this.multiSample.setProps(props.multiSample);
        if (props.drawPass) this.drawPass.setProps(props.drawPass);
    }

    get props(): ImageProps {
        return {
            transparentBackground: this._transparentBackground,
            postprocessing: this.postprocessing.props,
            multiSample: this.multiSample.props,
            drawPass: this.drawPass.props
        };
    }

    render() {
        Camera.copySnapshot(this._camera.state, this.camera.state);
        Viewport.set(this._camera.viewport, 0, 0, this._width, this._height);
        this._camera.update();

        this.renderer.setViewport(0, 0, this._width, this._height);

        if (this.multiSample.enabled) {
            this.multiSample.render(false, this._transparentBackground);
            this._colorTarget = this.multiSample.colorTarget;
        } else {
            this.drawPass.render(false, this._transparentBackground);
            if (this.postprocessing.enabled) {
                this.postprocessing.render(false);
                this._colorTarget = this.postprocessing.target;
            } else {
                this._colorTarget = this.drawPass.colorTarget;
            }
        }
    }

    getImageData(width: number, height: number) {
        this.setSize(width, height);
        this.render();
        const pd = this.colorTarget.getPixelData();
        return new ImageData(new Uint8ClampedArray(pd.array), pd.width, pd.height);
    }
}