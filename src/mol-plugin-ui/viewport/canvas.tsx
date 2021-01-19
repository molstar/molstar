/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';

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

    componentDidMount() {
        if (!this.canvas.current || !this.container.current || !this.plugin.initViewer(this.canvas.current!, this.container.current!)) {
            this.setState({ noWebGl: true });
            return;
        }
        this.handleLogo();
        this.subscribe(this.plugin.canvas3d!.reprCount, this.handleLogo);
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