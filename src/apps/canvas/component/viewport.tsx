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

interface ViewportProps {
    app: App
}

interface ViewportState {
    noWebGl: boolean
    pickingInfo: string
    taskInfo: string
}

export class Viewport extends React.Component<ViewportProps, ViewportState> {
    private container: HTMLDivElement | null = null;
    private canvas: HTMLCanvasElement | null = null;

    state: ViewportState = {
        noWebGl: false,
        pickingInfo: '',
        taskInfo: ''
    };

    handleResize() {
        this.props.app.canvas3d.handleResize()
    }

    componentDidMount() {
        if (!this.canvas || !this.container || !this.props.app.initViewer(this.canvas, this.container)) {
            this.setState({ noWebGl: true });
        }
        this.handleResize()

        const viewer = this.props.app.canvas3d

        viewer.input.resize.subscribe(() => this.handleResize())

        let prevLoci: Loci = EmptyLoci
        viewer.input.move.subscribe(async ({x, y, inside, buttons}) => {
            if (!inside || buttons) return
            const p = await viewer.identify(x, y)
            if (p) {
                const loci = viewer.getLoci(p)

                if (!areLociEqual(loci, prevLoci)) {
                    viewer.mark(prevLoci, MarkerAction.RemoveHighlight)
                    viewer.mark(loci, MarkerAction.Highlight)
                    prevLoci = loci

                    const label = labelFirst(loci)
                    const pickingInfo = `${label}`
                    this.setState({ pickingInfo })
                }
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