/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { resizeCanvas } from '../../mol-canvas3d/util';
import { Subject } from 'rxjs';
import { debounceTime } from 'rxjs/internal/operators/debounceTime';
import { PluginConfig } from '../../mol-plugin/config';
import { Color } from '../../mol-util/color';

interface ViewportCanvasState {
    noWebGl: boolean
    showLogo: boolean
}

export interface ViewportCanvasParams {
    logo?: React.FC,
    noWebGl?: React.FC,

    parentClassName?: string,
    parentStyle?: React.CSSProperties,
    hostClassName?: string,
    hostStyle?: React.CSSProperties,
}

export class ViewportCanvas extends PluginUIComponent<ViewportCanvasParams, ViewportCanvasState> {
    private container = React.createRef<HTMLDivElement>();
    private canvas = React.createRef<HTMLCanvasElement>();

    state: ViewportCanvasState = {
        noWebGl: false,
        showLogo: true
    };

    private handleLogo = () => {
        this.setState({ showLogo: !this.plugin.canvas3d?.reprCount.value });
    }

    private handleResize = () => {
        const container = this.container.current;
        const canvas = this.canvas.current;
        if (container && canvas) {
            const pixelScale = this.plugin.config.get(PluginConfig.General.PixelScale) || 1;
            resizeCanvas(canvas, container, pixelScale);
            const [r, g, b] = Color.toRgbNormalized(this.plugin.canvas3d!.props.renderer.backgroundColor);
            const a = this.plugin.canvas3d!.props.transparentBackground ? 0 : 1;
            this.plugin.canvas3d!.webgl.clear(r, g, b, a);
            this.plugin.canvas3d!.handleResize();
        }
    }

    componentDidMount() {
        if (!this.canvas.current || !this.container.current || !this.plugin.initViewer(this.canvas.current!, this.container.current!)) {
            this.setState({ noWebGl: true });
            return;
        }
        this.handleLogo();
        this.handleResize();

        const canvas3d = this.plugin.canvas3d!;
        this.subscribe(canvas3d.reprCount, this.handleLogo);

        const resized = new Subject();
        const resize = () => resized.next();

        this.subscribe(resized.pipe(debounceTime(1000 / 24)), () => this.handleResize());
        this.subscribe(canvas3d.input.resize, resize);
        this.subscribe(canvas3d.interaction.click, e => this.plugin.behaviors.interaction.click.next(e));
        this.subscribe(canvas3d.interaction.drag, e => this.plugin.behaviors.interaction.drag.next(e));
        this.subscribe(canvas3d.interaction.hover, e => this.plugin.behaviors.interaction.hover.next(e));
        this.subscribe(this.plugin.layout.events.updated, resize);
    }

    componentWillUnmount() {
        super.componentWillUnmount();
        // TODO viewer cleanup
    }

    renderMissing() {
        if (this.props.noWebGl) {
            const C = this.props.noWebGl;
            return <C />;
        }

        return <div className='msp-no-webgl'>
            <div>
                <p><b>WebGL does not seem to be available.</b></p>
                <p>This can be caused by an outdated browser, graphics card driver issue, or bad weather. Sometimes, just restarting the browser helps.</p>
                <p>For a list of supported browsers, refer to <a href='http://caniuse.com/#feat=webgl' target='_blank'>http://caniuse.com/#feat=webgl</a>.</p>
            </div>
        </div>;
    }

    render() {
        if (this.state.noWebGl) return this.renderMissing();

        const Logo = this.props.logo;

        return <div className={this.props.parentClassName || 'msp-viewport'} style={this.props.parentStyle}>
            <div className={this.props.hostClassName || 'msp-viewport-host3d'} style={this.props.hostStyle} ref={this.container}>
                <canvas ref={this.canvas} />
            </div>
            {(this.state.showLogo && Logo) && <Logo />}
        </div>;
    }
}