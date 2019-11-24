/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { CollapsableControls, CollapsableState } from './base';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ParameterControls } from './controls/parameters';
import { setCanvasSize } from '../../mol-canvas3d/util';

interface ImageControlsState extends CollapsableState {
    showPreview: boolean

    size: 'canvas' | 'custom'
    width: number
    height: number

    isDisabled: boolean
}

const maxWidthUi = 260
const maxHeightUi = 180

export class ImageControls<P, S extends ImageControlsState> extends CollapsableControls<P, S> {
    private canvasRef = React.createRef<HTMLCanvasElement>()

    private canvas: HTMLCanvasElement
    private canvasContext: CanvasRenderingContext2D

    // private imagePass: ImagePass

    get imagePass() {
        return this.plugin.helpers.viewportScreenshot!.imagePass;
    }

    constructor(props: P, context?: any) {
        super(props, context)

        this.subscribe(this.plugin.events.canvas3d.initialized, () => this.forceUpdate())
    }

    private getSize() {
        return this.state.size === 'canvas' ? {
            width: this.plugin.canvas3d.webgl.gl.drawingBufferWidth,
            height: this.plugin.canvas3d.webgl.gl.drawingBufferHeight
        } : {
            width: this.state.width,
            height: this.state.height
        }
    }

    private preview = () => {
        const { width, height } = this.getSize()
        if (width <= 0 || height <= 0) return

        let w: number, h: number
        const aH = maxHeightUi / height
        const aW = maxWidthUi / width
        if (aH < aW) {
            h = Math.round(Math.min(maxHeightUi, height))
            w = Math.round(width * (h / height))
        } else {
            w = Math.round(Math.min(maxWidthUi, width))
            h = Math.round(height * (w / width))
        }
        setCanvasSize(this.canvas, w, h)
        const { pixelRatio } = this.plugin.canvas3d.webgl
        const pw = Math.round(w * pixelRatio)
        const ph = Math.round(h * pixelRatio)
        const imageData = this.imagePass.getImageData(pw, ph)
        this.canvasContext.putImageData(imageData, 0, 0)
    }

    private download = () => {
        this.plugin.helpers.viewportScreenshot?.download();
    }

    private syncCanvas() {
        if (!this.canvasRef.current) return
        if (this.canvasRef.current === this.canvas) return

        this.canvas = this.canvasRef.current
        const ctx = this.canvas.getContext('2d')
        if (!ctx) throw new Error('Could not get canvas 2d context')
        this.canvasContext = ctx
    }

    private handlePreview() {
        if (this.state.showPreview) {
            this.syncCanvas()
            this.preview()
        }
    }

    componentDidUpdate() {
        this.plugin.helpers.viewportScreenshot!.size = this.getSize();
        this.handlePreview()
    }

    componentDidMount() {
        this.handlePreview()

        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => {
            this.imagePass.setProps({
                multiSample: { mode: 'on', sampleLevel: 2 },
                postprocessing: this.plugin.canvas3d.props.postprocessing
            })
            this.handlePreview()
        })

        this.subscribe(this.plugin.canvas3d.didDraw, () => {
            this.handlePreview()
        })

        this.subscribe(this.plugin.state.dataState.events.isUpdating, v => this.setState({ isDisabled: v }))
    }

    private togglePreview = () => this.setState({ showPreview: !this.state.showPreview })

    private setProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        if (p.name === 'size') {
            if (p.value.name === 'custom') {
                this.setState({ size: p.value.name, width: p.value.params.width, height: p.value.params.height })
            } else {
                this.setState({ size: p.value.name })
            }
        }
    }

    private get params () {
        const max = Math.min(this.plugin.canvas3d ? this.plugin.canvas3d.webgl.maxRenderbufferSize : 4096, 8192)
        const { width, height } = this.defaultState()
        return {
            size: PD.MappedStatic('custom', {
                canvas: PD.Group({}),
                custom: PD.Group({
                    width: PD.Numeric(width, { min: 128, max, step: 1 }),
                    height: PD.Numeric(height, { min: 128, max, step: 1 }),
                }, { isFlat: true })
            }, { options: [['canvas', 'Canvas'], ['custom', 'Custom']] })
        }
    }

    private get values () {
        return this.state.size === 'canvas'
            ? { size: { name: 'canvas', params: {} } }
            : { size: { name: 'custom', params: { width: this.state.width, height: this.state.height } } }
    }

    protected defaultState() {
        return {
            isCollapsed: false,
            header: 'Create Image',

            showPreview: false,

            size: 'canvas',
            width: 1920,
            height: 1080,

            isDisabled: false
        } as S
    }

    protected renderControls() {
        return <div>
            <div className='msp-control-row'>
                <button className='msp-btn msp-btn-block' onClick={this.download} disabled={this.state.isDisabled}>Download</button>
            </div>
            <ParameterControls params={this.params} values={this.values} onChange={this.setProps} isDisabled={this.state.isDisabled} />
            <div className='msp-control-group-wrapper'>
                <div className='msp-control-group-header'>
                    <button className='msp-btn msp-btn-block' onClick={this.togglePreview}>
                        <span className={`msp-icon msp-icon-${this.state.showPreview ? 'collapse' : 'expand'}`} />
                        Preview
                    </button>
                </div>
                {this.state.showPreview && <div className='msp-control-offset'>
                    <div className='msp-image-preview'>
                        <canvas width='0px' height='0px' ref={this.canvasRef} />
                    </div>
                </div>}
            </div>
        </div>
    }
}