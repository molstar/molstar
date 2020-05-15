/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { Subject } from 'rxjs';
import { debounceTime } from 'rxjs/operators';
import { resizeCanvas } from '../mol-canvas3d/util';
import { PluginCommands } from '../mol-plugin/commands';
import { PluginConfig } from '../mol-plugin/config';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { PluginUIComponent } from './base';
import { ControlGroup, IconButton } from './controls/common';
import { ToggleSelectionModeButton } from './structure/selection';
import { DownloadScreenshotControls } from './viewport/screenshot';
import { SimpleSettingsControl } from './viewport/simple-settings';
import { AutorenewSvg, CameraOutlinedSvg, BuildOutlinedSvg, FullscreenSvg, TuneSvg, CloseSvg } from './controls/icons';

interface ViewportControlsState {
    isSettingsExpanded: boolean,
    isScreenshotExpanded: boolean
}

interface ViewportControlsProps {
}

export class ViewportControls extends PluginUIComponent<ViewportControlsProps, ViewportControlsState> {
    private allCollapsedState: ViewportControlsState = {
        isSettingsExpanded: false,
        isScreenshotExpanded: false
    };

    state = { ...this.allCollapsedState } as ViewportControlsState;

    resetCamera = () => {
        PluginCommands.Camera.Reset(this.plugin, {});
    }

    private toggle(panel: keyof ViewportControlsState) {
        return (e?: React.MouseEvent<HTMLButtonElement>) => {
            this.setState({ ...this.allCollapsedState, [panel]: !this.state[panel] });
            e?.currentTarget.blur();
        };
    }

    toggleSettingsExpanded = this.toggle('isSettingsExpanded');
    toggleScreenshotExpanded = this.toggle('isScreenshotExpanded');

    toggleControls = () => {
        PluginCommands.Layout.Update(this.plugin, { state: { showControls: !this.plugin.layout.state.showControls } });
    }

    toggleExpanded = () => {
        PluginCommands.Layout.Update(this.plugin, { state: { isExpanded: !this.plugin.layout.state.isExpanded } });
    }

    setSettings = (p: { param: PD.Base<any>, name: string, value: any }) => {
        PluginCommands.Canvas3D.SetSettings(this.plugin, { settings: { [p.name]: p.value } });
    }

    setLayout = (p: { param: PD.Base<any>, name: string, value: any }) => {
        PluginCommands.Layout.Update(this.plugin, { state: { [p.name]: p.value } });
    }

    screenshot = () => {
        this.plugin.helpers.viewportScreenshot?.download();
    }

    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
        this.subscribe(this.plugin.layout.events.updated, () => this.forceUpdate());
    }

    icon(icon: React.FC, onClick: (e: React.MouseEvent<HTMLButtonElement>) => void, title: string, isOn = true) {
        return <IconButton svg={icon} toggleState={isOn} onClick={onClick} title={title} style={{ background: 'transparent' }} />;
    }

    onMouseMove = (e: React.MouseEvent) => {
        // ignore mouse moves when no button is held
        if (e.buttons === 0) e.stopPropagation();
    }

    render() {
        return <div className={'msp-viewport-controls'} onMouseMove={this.onMouseMove}>
            <div className='msp-viewport-controls-buttons'>
                <div>
                    <div className='msp-semi-transparent-background' />
                    {this.icon(AutorenewSvg, this.resetCamera, 'Reset Camera')}
                </div>
                <div>
                    <div className='msp-semi-transparent-background' />
                    {this.icon(CameraOutlinedSvg, this.toggleScreenshotExpanded, 'Screenshot / State Snapshot', this.state.isScreenshotExpanded)}
                </div>
                <div>
                    <div className='msp-semi-transparent-background' />
                    {this.icon(BuildOutlinedSvg, this.toggleControls, 'Toggle Controls Panel', this.plugin.layout.state.showControls)}
                    {this.plugin.config.get(PluginConfig.Viewport.ShowExpand) && this.icon(FullscreenSvg, this.toggleExpanded, 'Toggle Expanded Viewport', this.plugin.layout.state.isExpanded)}
                    {this.icon(TuneSvg, this.toggleSettingsExpanded, 'Settings / Controls Info', this.state.isSettingsExpanded)}
                </div>
                {this.plugin.config.get(PluginConfig.Viewport.ShowSelectionMode) && <div>
                    <div className='msp-semi-transparent-background' />
                    <ToggleSelectionModeButton />
                </div>}
            </div>
            {this.state.isScreenshotExpanded && <div className='msp-viewport-controls-panel'>
                <ControlGroup header='Screenshot / State' title='Click to close.' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleScreenshotExpanded}
                    topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                    <DownloadScreenshotControls close={this.toggleScreenshotExpanded} />
                </ControlGroup>
            </div>}
            {this.state.isSettingsExpanded && <div className='msp-viewport-controls-panel'>
                <ControlGroup header='Settings / Controls Info' title='Click to close.' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleSettingsExpanded}
                    topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                    <SimpleSettingsControl />
                </ControlGroup>
            </div>}
        </div>;
    }
}

export const Logo = () =>
    <a className='msp-logo' href='https://molstar.org' target='_blank' />;

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
        this.setState({ showLogo: !this.plugin.canvas3d?.reprCount.value });
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

        const resized = new Subject();
        const resize = () => resized.next();

        this.subscribe(resized.pipe(debounceTime(1000 / 24)), () => this.handleResize());
        this.subscribe(canvas3d.input.resize, resize);
        this.subscribe(canvas3d.interaction.click, e => this.plugin.behaviors.interaction.click.next(e));
        this.subscribe(canvas3d.interaction.hover, e => this.plugin.behaviors.interaction.hover.next(e));
        this.subscribe(this.plugin.layout.events.updated, resize);
    }

    componentWillUnmount() {
        super.componentWillUnmount();
        // TODO viewer cleanup
    }

    renderMissing() {
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

        return <div className='msp-viewport'>
            <div className='msp-viewport-host3d' ref={this.container}>
                <canvas ref={this.canvas} />
            </div>
            {this.state.showLogo && <Logo />}
        </div>;
    }
}