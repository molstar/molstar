/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Chetan Mishra <chetan.s115@gmail.com>
 */

import * as React from 'react';
import { throttleTime } from 'rxjs';
import { PluginCommands } from '../mol-plugin/commands';
import { PluginConfig } from '../mol-plugin/config';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { PluginUIComponent } from './base';
import { Button, ControlGroup, IconButton } from './controls/common';
import { AspectRatioSvg, AutorenewSvg, BuildOutlinedSvg, CameraOutlinedSvg, CloseSvg, FullscreenSvg, HeadsetVRSvg, LightModeSvg, TuneSvg } from './controls/icons';
import { ToggleSelectionModeButton } from './structure/selection';
import { ViewportCanvas } from './viewport/canvas';
import { DownloadScreenshotControls } from './viewport/screenshot';
import { SimpleSettingsControl } from './viewport/simple-settings';

interface ViewportControlsState {
    isSettingsExpanded: boolean,
    isScreenshotExpanded: boolean,
    isCameraResetEnabled: boolean
}

interface ViewportControlsProps {
}

export class ViewportControls extends PluginUIComponent<ViewportControlsProps, ViewportControlsState> {
    private allCollapsedState = {
        isSettingsExpanded: false,
        isScreenshotExpanded: false,
    };

    state: ViewportControlsState = {
        ...this.allCollapsedState,
        isCameraResetEnabled: true,
    };

    resetCamera = () => {
        PluginCommands.Camera.Reset(this.plugin, {});
    };

    private toggle(panel: keyof ViewportControlsState) {
        return (e?: React.MouseEvent<HTMLButtonElement>) => {
            this.setState(old => ({ ...old, ...this.allCollapsedState, [panel]: !this.state[panel] }));
            e?.currentTarget.blur();
        };
    }

    toggleSettingsExpanded = this.toggle('isSettingsExpanded');
    toggleScreenshotExpanded = this.toggle('isScreenshotExpanded');

    toggleControls = () => {
        PluginCommands.Layout.Update(this.plugin, { state: { showControls: !this.plugin.layout.state.showControls } });
    };

    toggleExpanded = () => {
        PluginCommands.Layout.Update(this.plugin, {
            state: {
                isExpanded: !this.plugin.layout.state.isExpanded,
                expandToFullscreen: false
            }
        });
    };

    toggleFullscreen = () => {
        PluginCommands.Layout.Update(this.plugin, {
            state: {
                expandToFullscreen: !this.plugin.layout.state.expandToFullscreen,
            }
        });
    };

    toggleXR = () => {
        if (this.plugin.canvas3d) {
            if (this.plugin.canvas3d.xr.isPresenting.value) {
                this.plugin.canvas3d.xr.end();
            } else {
                this.plugin.canvas3d.xr.request();
            }
        }
    };

    toggleIllumination = () => {
        if (!this.plugin.canvas3d) return;

        PluginCommands.Canvas3D.SetSettings(this.plugin, {
            settings: {
                illumination: {
                    ...this.plugin.canvas3d.props.illumination,
                    enabled: !this.plugin.canvas3d.props.illumination.enabled
                }
            }
        });
    };

    setSettings = (p: { param: PD.Base<any>, name: string, value: any }) => {
        PluginCommands.Canvas3D.SetSettings(this.plugin, { settings: { [p.name]: p.value } });
    };

    setLayout = (p: { param: PD.Base<any>, name: string, value: any }) => {
        PluginCommands.Layout.Update(this.plugin, { state: { [p.name]: p.value } });
    };

    screenshot = () => {
        this.plugin.helpers.viewportScreenshot?.download();
    };

    enableCameraReset = (enable: boolean) => {
        this.setState(old => ({ ...old, isCameraResetEnabled: enable }));
    };

    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
        this.subscribe(this.plugin.layout.events.updated, () => this.forceUpdate());
        if (this.plugin.canvas3d) {
            this.subscribe(
                this.plugin.canvas3d.camera.stateChanged.pipe(throttleTime(500, undefined, { leading: true, trailing: true })),
                snapshot => this.enableCameraReset(snapshot.radius !== 0 && snapshot.radiusMax !== 0)
            );
            this.subscribe(this.plugin.canvas3d.xr.isSupported, () => this.forceUpdate());
            this.subscribe(this.plugin.canvas3d.xr.isPresenting, () => this.forceUpdate());
        }
    }

    icon(icon: React.FC, onClick: (e: React.MouseEvent<HTMLButtonElement>) => void, title: string, isOn = true, disabled = false) {
        return <IconButton svg={icon} toggleState={isOn} onClick={onClick} title={title} style={{ background: 'transparent' }} disabled={disabled} />;
    }

    render() {
        const showXr = this.plugin.config.get(PluginConfig.Viewport.ShowXR);
        const xrIsSupported = !!this.plugin.canvas3d?.xr.isSupported.value;
        const xrIsPresenting = !!this.plugin.canvas3d?.xr.isPresenting.value;
        const xr = showXr === 'always' || (showXr === 'auto' && xrIsSupported);
        const xrTitle = !xrIsSupported ? 'Augmented/Virtual Reality unavailable' : (xrIsPresenting ? 'Exit XR' : 'Enter XR');
        const layoutState = this.plugin.layout.state;

        return <div className={'msp-viewport-controls'}>
            <div className='msp-viewport-controls-buttons'>
                {this.plugin.config.get(PluginConfig.Viewport.ShowReset) &&
                    <div className='msp-hover-box-wrapper'>
                        <div className='msp-semi-transparent-background' />
                        {this.icon(AutorenewSvg, this.resetCamera, 'Reset Zoom')}
                        <div className='msp-hover-box-body'>
                            <div className='msp-flex-column'>
                                <div className='msp-flex-row'>
                                    <Button onClick={() => this.resetCamera()} disabled={!this.state.isCameraResetEnabled} title='Set camera zoom to fit the visible scene into view'>
                                        Reset Zoom
                                    </Button>
                                </div>
                                <div className='msp-flex-row'>
                                    <Button onClick={() => PluginCommands.Camera.OrientAxes(this.plugin)} disabled={!this.state.isCameraResetEnabled} title='Align principal component axes of the loaded structures to the screen axes (“lay flat”)'>
                                        Orient Axes
                                    </Button>
                                </div>
                                <div className='msp-flex-row'>
                                    <Button onClick={() => PluginCommands.Camera.ResetAxes(this.plugin)} disabled={!this.state.isCameraResetEnabled} title='Align Cartesian axes to the screen axes'>
                                        Reset Axes
                                    </Button>
                                </div>
                            </div>
                        </div>
                        <div className='msp-hover-box-spacer'></div>
                    </div>
                }
                {this.plugin.config.get(PluginConfig.Viewport.ShowScreenshotControls) && <div>
                    <div className='msp-semi-transparent-background' />
                    {this.icon(CameraOutlinedSvg, this.toggleScreenshotExpanded, 'Screenshot / State Snapshot', this.state.isScreenshotExpanded)}
                </div>}
                <div>
                    <div className='msp-semi-transparent-background' />
                    {this.plugin.config.get(PluginConfig.Viewport.ShowControls) && this.icon(BuildOutlinedSvg, this.toggleControls, 'Toggle Controls Panel', this.plugin.layout.state.showControls)}
                    {this.plugin.config.get(PluginConfig.Viewport.ShowExpand) && <div className='msp-hover-box-wrapper'>
                        {this.icon(FullscreenSvg, this.toggleExpanded, 'Toggle Expanded Viewport', this.plugin.layout.state.isExpanded)}
                        <div className='msp-hover-box-body'>
                            <div className='msp-flex-column'>
                                <div className='msp-flex-row'>
                                    <Button onClick={this.toggleFullscreen}>
                                        {layoutState.expandToFullscreen ? 'Exit Fullscreen' : 'Fullscreen'}
                                    </Button>
                                </div>
                            </div>
                        </div>
                        <div className='msp-hover-box-spacer'></div>
                    </div>}
                    {!this.plugin.config.get(PluginConfig.Viewport.ShowExpand)
                        && this.plugin.config.get(PluginConfig.Viewport.ShowToggleFullscreen)
                        && this.icon(AspectRatioSvg, this.toggleFullscreen, 'Toggle Full Screen', this.plugin.layout.state.expandToFullscreen)}
                    {this.plugin.config.get(PluginConfig.Viewport.ShowSettings) && this.icon(TuneSvg, this.toggleSettingsExpanded, 'Settings / Controls Info', this.state.isSettingsExpanded)}
                    {this.plugin.config.get(PluginConfig.Viewport.ShowIllumination) && this.icon(LightModeSvg, this.toggleIllumination, 'Illumination', this.plugin.canvas3d?.props.illumination.enabled || false)}
                    {xr && this.icon(HeadsetVRSvg, this.toggleXR, xrTitle, xrIsPresenting, !xrIsSupported)}
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

export const Viewport = () => <ViewportCanvas logo={Logo} />;