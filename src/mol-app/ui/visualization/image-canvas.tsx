/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'

type State = { imageData: ImageData, width: number, height: number }

function getExtend(aspectRatio: number, maxWidth: number, maxHeight: number) {
    let width = maxWidth
    let height = width / aspectRatio
    if (height > maxHeight) {
        height = maxHeight
        width = height * aspectRatio
    }
    return { width, height }
}

export class ImageCanvas extends React.Component<{ imageData: ImageData, aspectRatio: number, maxWidth: number, maxHeight: number }, State> {
    private canvas: HTMLCanvasElement | null = null;
    private ctx: CanvasRenderingContext2D | null = null;

    updateStateFromProps() {
        this.setState({
            imageData: this.props.imageData,
            ...getExtend(this.props.aspectRatio, this.props.maxWidth, this.props.maxHeight)
        })
    }

    updateImage() {
        if (this.canvas) {
            this.canvas.width = this.state.imageData.width
            this.canvas.height = this.state.imageData.height
        }
        if (this.ctx) {
            this.ctx.putImageData(this.state.imageData, 0, 0)
        }
    }

    componentWillMount() {
        this.updateStateFromProps()
    }

    componentDidMount() {
        if (this.canvas && !this.ctx) {
            this.ctx = this.canvas.getContext('2d')
            if (this.ctx) this.ctx.imageSmoothingEnabled = false
        }
        this.updateImage()
    }

    componentWillReceiveProps() {
        this.updateStateFromProps()
    }

    componentDidUpdate() {
        this.updateImage()
    }

    render() {
        return <div
            className='molstar-image-canvas'
            style={{
                width: this.state.width + 6,
                height: this.state.height + 6,
                position: 'absolute',
                border: '3px white solid',
                bottom: 10,
                left: 10,
            }}
        >
            <canvas
                ref={elm => this.canvas = elm}
                style={{
                    width: this.state.width,
                    height: this.state.height,
                    imageRendering: 'pixelated'
                }}
            />
        </div>;
    }
}