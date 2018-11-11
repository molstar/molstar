/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { ButtonsType } from 'mol-util/input/input-observer';
import { Canvas3dIdentifyHelper } from 'mol-plugin/util/canvas3d-identify';
import { PluginComponent } from './base';

interface ViewportState {
    noWebGl: boolean
}

export class Viewport extends PluginComponent<{ }, ViewportState> {
    private container: HTMLDivElement | null = null;
    private canvas: HTMLCanvasElement | null = null;

    state: ViewportState = {
        noWebGl: false
    };

    private handleResize = () => {
         this.context.canvas3d.handleResize();
    }

    componentDidMount() {
        if (!this.canvas || !this.container || !this.context.initViewer(this.canvas, this.container)) {
            this.setState({ noWebGl: true });
        }
        this.handleResize();

        const canvas3d = this.context.canvas3d;
        this.subscribe(canvas3d.input.resize, this.handleResize);

        const idHelper = new Canvas3dIdentifyHelper(this.context, 15);

        this.subscribe(canvas3d.input.move, ({x, y, inside, buttons}) => {
            if (!inside || buttons) { return; }
            idHelper.move(x, y);
        });

        this.subscribe(canvas3d.input.leave, () => {
            idHelper.leave();
        });

        this.subscribe(canvas3d.input.click, ({x, y, buttons}) => {
            if (buttons !== ButtonsType.Flag.Primary) return;
            idHelper.select(x, y);
        });
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
        </div>;
    }
}