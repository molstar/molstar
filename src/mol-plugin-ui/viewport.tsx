/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Autorenew, BuildOutlined, CameraAltOutlined, Close, Crop, Fullscreen, Tune, CameraOutlined } from '@material-ui/icons';
import * as React from 'react';
import { resizeCanvas } from '../mol-canvas3d/util';
import { PluginCommands } from '../mol-plugin/commands';
import { PluginConfig } from '../mol-plugin/config';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { PluginUIComponent } from './base';
import { ControlGroup, IconButton } from './controls/common';
import { DownloadScreenshotControls } from './viewport/screenshot';
import { SimpleSettingsControl } from './viewport/simple-settings';

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

    toggleSelectionMode = () => {
        this.plugin.selectionMode = !this.plugin.selectionMode;
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
        this.subscribe(this.plugin.behaviors.interaction.selectionMode, () => this.forceUpdate());
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
                    {this.icon(Autorenew, this.resetCamera, 'Reset Camera')}
                </div>
                <div>
                    <div className='msp-semi-transparent-background' />
                    {this.icon(CameraOutlined, this.toggleScreenshotExpanded, 'Screenshot / State Snapshot', this.state.isScreenshotExpanded)}
                </div>
                <div>
                    <div className='msp-semi-transparent-background' />
                    {this.icon(BuildOutlined, this.toggleControls, 'Toggle Controls', this.plugin.layout.state.showControls)}
                    {this.plugin.config.get(PluginConfig.Viewport.ShowExpand) && this.icon(Fullscreen, this.toggleExpanded, 'Toggle Expanded', this.plugin.layout.state.isExpanded)}
                    {this.icon(Tune, this.toggleSettingsExpanded, 'Settings / Controls Info', this.state.isSettingsExpanded)}
                </div>
                {this.plugin.config.get(PluginConfig.Viewport.ShowSelectionMode) && <div>
                    <div className='msp-semi-transparent-background' />
                    {this.icon(Crop, this.toggleSelectionMode, 'Toggle Selection Mode', this.plugin.behaviors.interaction.selectionMode.value)}
                </div>}
            </div>
            {this.state.isScreenshotExpanded && <div className='msp-viewport-controls-panel'>
                <ControlGroup header='Screenshot / State Snapshot' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleScreenshotExpanded}
                    topRightIcon={Close} noTopMargin>
                    <DownloadScreenshotControls close={this.toggleScreenshotExpanded} />
                </ControlGroup>
            </div>}
            {this.state.isSettingsExpanded && <div className='msp-viewport-controls-panel'>
                <ControlGroup header='Settings / Controls Info' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleSettingsExpanded}
                    topRightIcon={Close} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                    <SimpleSettingsControl />
                </ControlGroup>
            </div>}
        </div>;
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
    </div>;

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