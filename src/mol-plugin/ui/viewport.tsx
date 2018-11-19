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
import { PluginCommands } from 'mol-plugin/command';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ParameterControls } from './controls/parameters';
import { Canvas3DParams } from 'mol-canvas3d/canvas3d';

interface ViewportState {
    noWebGl: boolean
}

export class ViewportControls extends PluginComponent {
    state = {
        isSettingsExpanded: false,
        settings: PD.getDefaultValues(Canvas3DParams)
    }

    resetCamera = () => {
        PluginCommands.Camera.Reset.dispatch(this.plugin, {});
    }

    toggleSettingsExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isSettingsExpanded: !this.state.isSettingsExpanded });
        e.currentTarget.blur();
    }

    setSettings = (p: { param: PD.Base<any>, name: string, value: any }) => {
        this.plugin.canvas3d.setProps({ [p.name]: p.value })
        this.setState({ settings: this.plugin.canvas3d.props })
    }

    render() {
        return <div className={'msp-viewport-controls'}>
            <div className={'msp-row'}>
                <button className='msp-btn msp-btn-link' onClick={this.resetCamera}>↻ Camera</button>
                <div>
                    <button className='msp-btn msp-btn-link'onClick={this.toggleSettingsExpanded}>Settings</button>
                    <div className='msp-control-offset' style={{ display: this.state.isSettingsExpanded ? 'block' : 'none' }}>
                        <ParameterControls params={Canvas3DParams} values={this.state.settings} onChange={this.setSettings} />
                    </div>
                </div>
            </div>
        </div>
    }
}

export class Viewport extends PluginComponent<{ }, ViewportState> {
    private container: HTMLDivElement | null = null;
    private canvas: HTMLCanvasElement | null = null;

    state: ViewportState = {
        noWebGl: false
    };

    private handleResize = () => {
         this.plugin.canvas3d.handleResize();
    }

    componentDidMount() {
        if (!this.canvas || !this.container || !this.plugin.initViewer(this.canvas, this.container)) {
            this.setState({ noWebGl: true });
        }
        this.handleResize();

        const canvas3d = this.plugin.canvas3d;
        this.subscribe(canvas3d.input.resize, this.handleResize);

        const idHelper = new Canvas3dIdentifyHelper(this.plugin, 15);

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

        return <div className='msp-viewport'>
            <div className='msp-viewport-host3d' ref={elm => this.container = elm}>
                <canvas ref={elm => this.canvas = elm}></canvas>
            </div>
        </div>;
    }
}