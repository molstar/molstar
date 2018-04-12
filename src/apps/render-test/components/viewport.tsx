/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import State from '../state'

export default class Viewport extends React.Component<{ state: State }, { initialized: boolean }> {
    private canvasContainer: HTMLDivElement | null = null;
    state = { initialized: false }

    componentDidMount() {
        if (this.canvasContainer) this.props.state.initRenderer(this.canvasContainer).then(() => this.setState({ initialized: true }))
    }

    render() {
        return <div ref={elm => this.canvasContainer = elm} style={{ height: '100%' }}>
        </div>
    }
}