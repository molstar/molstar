/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginUIComponent } from './base';
import { PluginCommands } from '../../mol-plugin/command';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ParameterControls } from './controls/parameters';
import { Canvas3DParams } from '../../mol-canvas3d/canvas3d';
import { PluginLayoutStateParams } from '../../mol-plugin/layout';
import { ControlGroup, IconButton } from './controls/common';
import { resizeCanvas } from '../../mol-canvas3d/util';
import { Interactivity } from '../util/interactivity';

interface ViewportControlsState {
    isSettingsExpanded: boolean
}

interface ViewportControlsProps {
    hideSettingsIcon?: boolean
}

export class ViewportControls extends PluginUIComponent<ViewportControlsProps, ViewportControlsState> {
    state = {
        isSettingsExpanded: false,
    }

    resetCamera = () => {
        PluginCommands.Camera.Reset.dispatch(this.plugin, {});
    }

    toggleSettingsExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isSettingsExpanded: !this.state.isSettingsExpanded });
        e.currentTarget.blur();
    }

    toggleControls = () => {
        PluginCommands.Layout.Update.dispatch(this.plugin, { state: { showControls: !this.plugin.layout.state.showControls } });
    }

    toggleExpanded = () => {
        PluginCommands.Layout.Update.dispatch(this.plugin, { state: { isExpanded: !this.plugin.layout.state.isExpanded } });
    }

    setSettings = (p: { param: PD.Base<any>, name: string, value: any }) => {
        PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: { [p.name]: p.value } });
    }

    setLayout = (p: { param: PD.Base<any>, name: string, value: any }) => {
        PluginCommands.Layout.Update.dispatch(this.plugin, { state: { [p.name]: p.value } });
    }

    setInteractivityProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        PluginCommands.Interactivity.SetProps.dispatch(this.plugin, { props: { [p.name]: p.value } });
    }

    screenshot = () => {
        this.plugin.helpers.viewportScreenshot?.download();
    }

    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
        this.subscribe(this.plugin.layout.events.updated, () => this.forceUpdate());
        this.subscribe(this.plugin.events.interactivity.propsUpdated, () => this.forceUpdate());
    }

    icon(name: string, onClick: (e: React.MouseEvent<HTMLButtonElement>) => void, title: string, isOn = true) {
        return <IconButton icon={name} toggleState={isOn} onClick={onClick} title={title} />;
    }

    onMouseMove = (e: React.MouseEvent) => {
        // ignore mouse moves when no button is held
        if (e.buttons === 0) e.stopPropagation()
    }

    render() {
        return <div className={'msp-viewport-controls'} onMouseMove={this.onMouseMove}>
            <div className='msp-viewport-controls-buttons'>
                <div className='msp-semi-transparent-background' />
                {this.icon('reset-scene', this.resetCamera, 'Reset Camera')}
                {this.icon('screenshot', this.screenshot, 'Download Screenshot')}
                {this.icon('tools', this.toggleControls, 'Toggle Controls', this.plugin.layout.state.showControls)}
                {this.icon('expand-layout', this.toggleExpanded, 'Toggle Expanded', this.plugin.layout.state.isExpanded)}
                {!this.props.hideSettingsIcon && this.icon('settings', this.toggleSettingsExpanded, 'Settings', this.state.isSettingsExpanded)}
            </div>
            {this.state.isSettingsExpanded && <div className='msp-viewport-controls-scene-options'>
                <ControlGroup header='Layout' initialExpanded={true}>
                    <ParameterControls params={PluginLayoutStateParams} values={this.plugin.layout.state} onChange={this.setLayout} />
                </ControlGroup>
                <ControlGroup header='Interactivity' initialExpanded={true}>
                    <ParameterControls params={Interactivity.Params} values={this.plugin.interactivity.props} onChange={this.setInteractivityProps} />
                </ControlGroup>
                {this.plugin.canvas3d && <ControlGroup header='Viewport' initialExpanded={true}>
                    <ParameterControls params={Canvas3DParams} values={this.plugin.canvas3d.props} onChange={this.setSettings} />
                </ControlGroup>}
            </div>}
        </div>
    }
}

export const Logo = () =>
    <div className='msp-logo'>
        <div>
            <div>
                <div />
                <div className='msp-logo-image' />
            </div>
        </div>
    </div>

interface ViewportState {
    noWebGl: boolean
    showLogo: boolean
}

export class Viewport extends PluginUIComponent<{ }, ViewportState> {
    private container = React.createRef<HTMLDivElement>();
    private canvas = React.createRef<HTMLCanvasElement>();

    state: ViewportState = {
        noWebGl: false,
        showLogo: true
    };

    private handleLogo = () => {
        this.setState({ showLogo: !this.plugin.canvas3d?.reprCount.value })
    }

    private handleResize = () => {
        const container = this.container.current;
        const canvas = this.canvas.current;
        if (container && canvas) {
            resizeCanvas(canvas, container);
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
        this.subscribe(canvas3d.input.resize, this.handleResize);

        this.subscribe(canvas3d.interaction.click, e => this.plugin.behaviors.interaction.click.next(e));
        this.subscribe(canvas3d.interaction.hover, e => this.plugin.behaviors.interaction.hover.next(e));
        this.subscribe(this.plugin.layout.events.updated, () => {
            setTimeout(this.handleResize, 50);
        });
    }

    componentWillUnmount() {
        if (super.componentWillUnmount) super.componentWillUnmount();
        // TODO viewer cleanup
    }

    renderMissing() {
        return <div className='msp-no-webgl'>
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
            <div className='msp-viewport-host3d' ref={this.container}>
                <canvas ref={this.canvas} />
            </div>
            {this.state.showLogo && <Logo />}
        </div>;
    }
}