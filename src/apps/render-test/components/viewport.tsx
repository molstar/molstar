/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import State from '../state'

export default class Viewport extends React.Component<{ state: State }, { initialized: boolean }> {
    private container: HTMLDivElement | null = null;
    private canvas: HTMLCanvasElement | null = null;
    state = { initialized: false }

    componentDidMount() {
        if (this.container && this.canvas) {
            this.props.state.initRenderer(this.canvas, this.container).then(() => {
                this.setState({ initialized: true })
            })
        }
    }

    render() {
        return <div ref={elm => this.container = elm} style={{ height: '100%' }}>
            <canvas ref={elm => this.canvas = elm}></canvas>
        </div>
    }
}