/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'

import { ViewportController } from '../../controller/visualization/viewport'
import { View } from '../view';
import { HelpBox, Toggle, Button } from '../controls/common'
import { Slider } from '../controls/slider'
import { ImageCanvas } from './image-canvas';
import { InteractivityEvents } from '../../event/basic';
import { labelFirst } from 'mol-view/label';

export class ViewportControls extends View<ViewportController, { showSceneOptions?: boolean, showHelp?: boolean }, {}> {
    state = { showSceneOptions: false, showHelp: false };

    private help() {
        return <div className='molstar-viewport-controls-scene-options molstar-control'>
            <HelpBox title='Rotate' content={<div><div>Left button</div><div>One finger touch</div></div>} />
            <HelpBox title='Zoom' content={<div><div>Right button</div><div>Pinch</div></div>} />
            <HelpBox title='Move' content={<div><div>Middle button</div><div>Two finger touch</div></div>} />
            <HelpBox title='Slab' content={<div><div>Mouse wheel</div><div>Three finger touch</div></div>} />
        </div>
    }

    render() {
        let state = this.controller.latestState;

        let options: any;

        let layoutController = this.controller.context.layout;
        let layoutState = layoutController.latestState;
        if (this.state.showSceneOptions) {
            options = <div className='molstar-viewport-controls-scene-options molstar-control'>
                <Toggle onChange={v => this.controller.setState({ enableFog: v })} value={state.enableFog!} label='Fog' />
                <Slider label='FOV' min={30} max={90} onChange={v => this.controller.setState({ cameraFOV: v }) } value={state.cameraFOV!} />
                <Slider label='Camera Speed' min={1} max={10} step={0.01} onChange={v => this.controller.setState({ cameraSpeed: v }) } value={state.cameraSpeed!} />
            </div>;
        } else if (this.state.showHelp) {
            options = this.help();
        }

        let controlsShown = !layoutState.hideControls;
        return <div className='molstar-viewport-controls' onMouseLeave={() => this.setState({ showSceneOptions: false, showHelp: false })}>
            <div className='molstar-viewport-controls-buttons'>
                <Button
                    style='link'
                    active={this.state.showHelp}
                    customClass={'molstar-btn-link-toggle-' + (this.state.showHelp ? 'on' : 'off')}
                    icon='help-circle'
                    onClick={(e) => this.setState({ showHelp: !this.state.showHelp, showSceneOptions: false }) } title='Controls Help' />
                <Button
                    style='link'
                    active={this.state.showSceneOptions}
                    customClass={'molstar-btn-link-toggle-' + (this.state.showSceneOptions ? 'on' : 'off')}
                    icon='settings'
                    onClick={(e) => this.setState({ showSceneOptions: !this.state.showSceneOptions, showHelp: false }) } title='Scene Options' />
                <Button
                    style='link'
                    icon='screenshot'
                    onClick={(e) => this.controller.context.stage.viewer.downloadScreenshot()}
                    title='Screenshot' />
                <Button   onClick={() => { layoutController.update({ hideControls: controlsShown }); this.forceUpdate(); } }
                    icon='tools' title={controlsShown ? 'Hide Controls' : 'Show Controls'} active={controlsShown }
                    customClass={'molstar-btn-link-toggle-' + (controlsShown  ? 'on' : 'off')}
                    style='link' />
                <Button   onClick={() => layoutController.update({ isExpanded: !layoutState.isExpanded  }) }
                    icon='expand-layout' title={layoutState.isExpanded ? 'Collapse' : 'Expand'} active={layoutState.isExpanded }
                    customClass={'molstar-btn-link-toggle-' + (layoutState.isExpanded  ? 'on' : 'off')}
                    style='link' />
                <Button
                    style='link'
                    icon='reset-scene'
                    onClick={(e) => this.controller.context.stage.viewer.resetCamera()}
                    title='Reset camera' />
            </div>
            {options}
        </div>;
    }
}

export const Logo = () =>
    <div className='molstar-logo'>
        <div>
            <div>
                <div />
                <div className='molstar-logo-image' />
            </div>
        </div>
    </div>


type ViewportState = {
    noWebGl: boolean,
    showLogo: boolean,
    aspectRatio: number,
    width: number
    height: number
    images: { [k: string]: ImageData }
    info: string
}

export class Viewport extends View<ViewportController, ViewportState, { noWebGl?: boolean, showLogo?: boolean, aspectRatio: number, info: string }> {
    private container: HTMLDivElement | null = null;
    private canvas: HTMLCanvasElement | null = null;
    private defaultBg = { r: 1, g: 1, b: 1 }
    state: ViewportState = {
        noWebGl: false,
        showLogo: true,
        images: {},
        aspectRatio: 1,
        width: 0,
        height: 0,
        info: ''
    };

    handleResize() {
        if (this.container) {
            this.setState({
                aspectRatio: this.container.clientWidth / this.container.clientHeight,
                width: this.container.clientWidth,
                height: this.container.clientHeight
            })
        }
    }

    componentDidMount() {
        if (!this.canvas || !this.container || !this.controller.context.initStage(this.canvas, this.container)) {
            this.setState({ noWebGl: true });
        }
        this.handleResize()

        const viewer = this.controller.context.stage.viewer

        viewer.reprCount.subscribe(count => {
            this.setState({
                showLogo: false
                // showLogo: count === 0
            })
        })

        viewer.didDraw.subscribe(() => {
            // this.setState({ imageData: viewer.getImageData() })
            this.setState({
                images: {
                    'object': viewer.getImageData('pickObject'),
                    'instance': viewer.getImageData('pickInstance'),
                    'element': viewer.getImageData('pickElement')
                }
            })
        })

        viewer.input.resize.subscribe(() => this.handleResize())

        viewer.input.move.subscribe(({x, y, inside}) => {
            if (!inside) return
            const p = viewer.identify(x, y)
            const loci = viewer.getLoci(p)
            InteractivityEvents.HighlightLoci.dispatch(this.controller.context, loci);
            
            // TODO use LabelLoci event and make configurable
            const label = labelFirst(loci)
            const info = `Object: ${p.objectId}, Instance: ${p.instanceId}, Element: ${p.elementId}, Label: ${label}`
            this.setState({ info })
        })

        // TODO filter only for left button?
        viewer.input.click.subscribe(({x, y}) => {
            const loci = viewer.getLoci(viewer.identify(x, y))
            InteractivityEvents.SelectLoci.dispatch(this.controller.context, loci);
        })
    }

    componentWillUnmount() {
        super.componentWillUnmount();
        this.controller.context.destroy();
    }

    renderMissing() {
        return <div className='molstar-no-webgl'>
            <div>
                <p><b>WebGL does not seem to be available.</b></p>
                <p>This can be caused by an outdated browser, graphics card driver issue, or bad weather. Sometimes, just restarting the browser helps.</p>
                <p>For a list of supported browsers, refer to <a href='http://caniuse.com/#feat=webgl' target='_blank'>http://caniuse.com/#feat=webgl</a>.</p>
            </div>
        </div>
    }

    render() {
        if (this.state.noWebGl) return this.renderMissing();

        const color = this.controller.latestState.clearColor! || this.defaultBg;
        return <div className='molstar-viewport' style={{ backgroundColor: `rgb(${255 * color.r}, ${255 * color.g}, ${255 * color.b})` }}>
            <div ref={elm => this.container = elm} className='molstar-viewport-container'>
                <canvas ref={elm => this.canvas = elm} className='molstar-viewport-canvas'></canvas>
            </div>
            {this.state.showLogo ? <Logo /> : void 0}
            <ViewportControls controller={this.controller} />
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
                {this.state.info}
            </div>
            <div
                style={{
                    position: 'absolute',
                    bottom: 10,
                    left: 10,
                }}
            >
                {Object.keys(this.state.images).map(k => {
                    const imageData = this.state.images[k]
                    return <ImageCanvas
                        key={k}
                        imageData={imageData}
                        aspectRatio={this.state.aspectRatio}
                        maxWidth={this.state.width / 4}
                        maxHeight={this.state.height / 4}
                    />
                })}
            </div>
        </div>;
    }
}