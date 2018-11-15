/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { App } from '../app';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { EmptyLoci, Loci, areLociEqual } from 'mol-model/loci';
import { labelFirst } from 'mol-theme/label';
import { ButtonsType } from 'mol-util/input/input-observer';
import { throttleTime } from 'rxjs/operators'
import { Camera } from 'mol-canvas3d/camera';
import { ColorParamComponent } from 'mol-app/component/parameter/color';
import { Color } from 'mol-util/color';
import { ParamDefinition as PD } from 'mol-util/param-definition'

interface ViewportProps {
    app: App
}

interface ViewportState {
    noWebGl: boolean
    pickingInfo: string
    taskInfo: string
    cameraMode: Camera.Mode
    backgroundColor: Color
}

const BackgroundColorParam = PD.Color(Color(0x000000), { label: 'Background Color' })

export class Viewport extends React.Component<ViewportProps, ViewportState> {
    private container: HTMLDivElement | null = null;
    private canvas: HTMLCanvasElement | null = null;

    state: ViewportState = {
        noWebGl: false,
        pickingInfo: '',
        taskInfo: '',
        cameraMode: 'perspective',
        backgroundColor: Color(0x000000)
    };

    handleResize() {
        this.props.app.canvas3d.handleResize()
    }

    componentDidMount() {
        if (!this.canvas || !this.container || !this.props.app.initViewer(this.canvas, this.container)) {
            this.setState({ noWebGl: true });
        }
        this.handleResize()

        this.setState({
            cameraMode: this.props.app.canvas3d.props.cameraMode,
            backgroundColor: this.props.app.canvas3d.props.backgroundColor
        })

        const canvas3d = this.props.app.canvas3d

        canvas3d.input.resize.subscribe(() => this.handleResize())

        let prevHighlightLoci: Loci = EmptyLoci
        // TODO can the 'only ever have one extra element in the queue' functionality be done with rxjs?
        let highlightQueueLength = 0
        canvas3d.input.move.pipe(throttleTime(50)).subscribe(async ({x, y, inside, buttons}) => {
            if (!inside || buttons || highlightQueueLength > 2) return
            ++highlightQueueLength
            const p = await canvas3d.identify(x, y)
            --highlightQueueLength
            if (p) {
                const { loci } = canvas3d.getLoci(p)

                if (!areLociEqual(loci, prevHighlightLoci)) {
                    canvas3d.mark(prevHighlightLoci, MarkerAction.RemoveHighlight)
                    canvas3d.mark(loci, MarkerAction.Highlight)
                    prevHighlightLoci = loci

                    const label = labelFirst(loci)
                    const pickingInfo = `${label}`
                    this.setState({ pickingInfo })
                }
            }
        })

        canvas3d.input.click.subscribe(async ({x, y, buttons}) => {
            if (buttons !== ButtonsType.Flag.Primary) return
            const p = await canvas3d.identify(x, y)
            if (p) {
                const { loci } = canvas3d.getLoci(p)
                canvas3d.mark(loci, MarkerAction.Toggle)
            }
        })

        this.props.app.taskCountChanged.subscribe(({ count, info }) => {
            this.setState({ taskInfo: count > 0 ? info : '' })
        })
    }

    componentWillUnmount() {
        if (super.componentWillUnmount) super.componentWillUnmount();
        // TODO viewer cleanup
    }

    renderMissing() {
        return <div>
            <div>
                <p><b>WebGL does not seem to be available.</b></p>
                <p>This can be caused by an outdated browser, graphics card driver issue, or bad weather. Sometimes, just restarting the browser helps.</p>
                <p>For a list of supported browsers, refer to <a href='http://caniuse.com/#feat=webgl' target='_blank'>http://caniuse.com/#feat=webgl</a>.</p>
            </div>
        </div>
    }

    render() {
        if (this.state.noWebGl) return this.renderMissing();

        return <div style={{ backgroundColor: 'rgb(0, 0, 0)', width: '100%', height: '100%'}}>
            <div ref={elm => this.container = elm} style={{width: '100%', height: '100%'}}>
                <canvas ref={elm => this.canvas = elm}></canvas>
            </div>
            <div
                style={{
                    position: 'absolute',
                    top: 10,
                    left: 10,
                    padding: 10,
                    color: 'lightgrey',
                    background: 'rgba(0, 0, 0, 0.2)'
                }}
            >
                {this.state.pickingInfo}
            </div>
            <div
                style={{
                    position: 'absolute',
                    bottom: 10,
                    right: 10,
                    padding: 10,
                    color: 'lightgrey',
                    background: 'rgba(0, 0, 0, 0.2)'
                }}
            >
                <div>
                    <span>Camera Mode </span>
                    <select
                        value={this.state.cameraMode}
                        style={{width: '150'}}
                        onChange={e => {
                            const p = { cameraMode: e.target.value as Camera.Mode }
                            this.props.app.canvas3d.setProps(p)
                            this.setState(p)
                        }}
                    >
                        <option value='perspective'>Perspective</option>
                        <option value='orthographic'>Orthographic</option>
                    </select>
                </div>
                <ColorParamComponent
                    label={BackgroundColorParam.label || ''}
                    param={BackgroundColorParam}
                    value={this.state.backgroundColor}
                    onChange={value => {
                        const p = { backgroundColor: value }
                        this.props.app.canvas3d.setProps(p)
                        this.setState(p)
                    }}
                />
            </div>
            { this.state.taskInfo ?
                <div
                    style={{
                        position: 'absolute',
                        top: 10,
                        right: 10,
                        padding: 10,
                        color: 'lightgrey',
                        background: 'rgba(0, 0, 0, 0.2)'
                    }}
                >
                    {this.state.taskInfo}
                </div>
            : '' }
        </div>;
    }
}