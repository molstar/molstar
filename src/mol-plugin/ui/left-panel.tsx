/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginUIComponent } from './base';
import { StateTree } from './state/tree';
import { IconButton, SectionHeader, ControlGroup } from './controls/common';
import { StateObjectActions } from './state/actions';
import { StateTransform } from '../../mol-state';
import { PluginCommands } from '../command';
import { ParameterControls } from './controls/parameters';
import { Canvas3DParams } from '../../mol-canvas3d/canvas3d';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StateSnapshots } from './state/snapshots';

type TabName = 'root' | 'data' | 'behavior' | 'states' | 'viewport-settings'

export class LeftPanelControls extends PluginUIComponent<{}, { currentTab: TabName }> {
    state = { currentTab: 'data' as TabName };

    componentDidMount() {
        // this.subscribe(this.plugin.state.behavior.kind, () => this.forceUpdate());
    }

    set(kind: TabName) {
        switch (kind) {
            case 'data': this.plugin.state.setKind('data'); break;
            case 'behavior': this.plugin.state.setKind('behavior'); break;
        }
        this.setState({ currentTab: kind });
    }

    tabs: { [K in TabName]: JSX.Element } = {
        'root': <>
            <SectionHeader icon='home' title='Home' />
            <StateObjectActions state={this.plugin.state.dataState} nodeRef={StateTransform.RootRef} hideHeader={true} initiallyCollapsed={true} alwaysExpandFirst={true} />
        </>,
        'data': <>
        <SectionHeader icon='flow-tree' title='State Tree' />
            <StateTree state={this.plugin.state.dataState} />
        </>,
        'states': <StateSnapshots />,
        'behavior': <>
            <SectionHeader icon='address' title='Plugin Behavior' />
            <StateTree state={this.plugin.state.behaviorState} />
        </>,
        'viewport-settings': <>
            <SectionHeader icon='settings' title='Plugin Settings' />
            <FullSettings />
        </>
    }

    render() {
        const tab = this.state.currentTab;

        return <div className='msp-left-panel-controls'>
            <div className='msp-left-panel-controls-buttons'>
                <IconButton icon='home' toggleState={tab === 'root'} onClick={() => this.set('root')} title='Home' />
                <IconButton icon='flow-tree' toggleState={tab === 'data'} onClick={() => this.set('data')} title='State Tree' />
                <IconButton icon='floppy' toggleState={tab === 'states'} onClick={() => this.set('states')} title='Plugin State' />
                <div className='msp-left-panel-controls-buttons-bottom'>
                    <IconButton icon='address' toggleState={tab === 'behavior'} onClick={() => this.set('behavior')} title='Plugin Behavior' />
                    <IconButton icon='settings' toggleState={tab === 'viewport-settings'} onClick={() => this.set('viewport-settings')} title='Viewport Settings' />
                </div>
            </div>
            <div className='msp-scrollable-container'>
                {this.tabs[tab]}
            </div>
        </div>;
    }
}

class FullSettings extends PluginUIComponent {
    private setSettings = (p: { param: PD.Base<any>, name: string, value: any }) => {
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

    render() {
        return <>
            {/* <ControlGroup header='Layout' initialExpanded={true}>
                <ParameterControls params={PluginLayoutStateParams} values={this.plugin.layout.state} onChange={this.setLayout} />
            </ControlGroup>
            <ControlGroup header='Interactivity' initialExpanded={true}>
                <ParameterControls params={Interactivity.Params} values={this.plugin.interactivity.props} onChange={this.setInteractivityProps} />
            </ControlGroup> */}
            {this.plugin.canvas3d && <ControlGroup header='Viewport' initialExpanded={true}>
                <ParameterControls params={Canvas3DParams} values={this.plugin.canvas3d.props} onChange={this.setSettings} />
            </ControlGroup>}
        </>
    }
}
