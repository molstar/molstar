/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mp4EncoderUI } from '../../../extensions/mp4-export/ui';
import { CollapsableControls, CollapsableState, PluginUIComponent } from '../../../mol-plugin-ui/base';
import { SectionHeader } from '../../../mol-plugin-ui/controls/common';
import { ParameterControls } from '../../../mol-plugin-ui/controls/parameters';
import { PluginCommands } from '../../../mol-plugin/commands';
import { StructureMeasurementsControls } from '../../../mol-plugin-ui/structure/measurements';
import { MesoscaleExplorerState } from '../app';
import { MesoscaleState } from '../data/state';
import { EntityControls, FocusInfo, ModelInfo, SelectionInfo } from './entities';
import { LoaderControls, ExampleControls, SessionControls, SnapshotControls, DatabaseControls, MesoQuickStylesControls, ExplorerInfo } from './states';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { TuneSvg } from '../../../mol-plugin-ui/controls/icons';
import { RendererParams } from '../../../mol-gl/renderer';
import { TrackballControlsParams } from '../../../mol-canvas3d/controls/trackball';

const Spacer = () => <div style={{ height: '2em' }} />;

const ViewportParams = {
    renderer: PD.Group(RendererParams),
    trackball: PD.Group(TrackballControlsParams),
};

class ViewportSettingsUI extends CollapsableControls<{}, {}> {
    protected defaultState(): CollapsableState {
        return {
            header: 'Viewport Settings',
            isCollapsed: true,
            brand: { accent: 'cyan', svg: TuneSvg }
        };
    }

    protected renderControls(): JSX.Element | null {
        return <>
            {this.plugin.canvas3d && this.plugin.canvas3dContext && <>
                <ParameterControls params={ViewportParams} values={this.plugin.canvas3d.props} onChange={this.setSettings} />
            </>}
        </>;
    }

    private setSettings = (p: { param: PD.Base<any>, name: string, value: any }) => {
        PluginCommands.Canvas3D.SetSettings(this.plugin, { settings: { [p.name]: p.value } });
    };

    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
        this.subscribe(this.plugin.layout.events.updated, () => this.forceUpdate());
    }
}

export class LeftPanel extends PluginUIComponent {
    render() {
        const customState = this.plugin.customState as MesoscaleExplorerState;

        return <div className='msp-scrollable-container'>
            {customState.driver && <>
                <ExplorerInfo />
                <Spacer />
            </>}
            <SectionHeader title='Database' />
            <DatabaseControls />
            <Spacer />

            <SectionHeader title='Open' />
            <LoaderControls />
            <Spacer />

            {customState.examples?.length && <>
                <SectionHeader title='Example' />
                <ExampleControls />
                <Spacer />
            </>}

            <SectionHeader title='Session' />
            <SessionControls />
            <Spacer />

            <SectionHeader title='Snapshots' />
            <SnapshotControls />
            <Spacer />

            <Mp4EncoderUI />
            <ViewportSettingsUI />
        </div>;
    }
}

export class RightPanel extends PluginUIComponent<{}, { isDisabled: boolean }> {
    state = {
        isDisabled: false,
    };

    get hasModelInfo() {
        return (
            MesoscaleState.has(this.plugin) &&
            !!(MesoscaleState.get(this.plugin).description ||
                MesoscaleState.get(this.plugin).link)
        );
    }

    get hasFocusInfo() {
        return (
            MesoscaleState.has(this.plugin) &&
            !!(MesoscaleState.get(this.plugin).focusInfo !== '')
        );
    }

    componentDidMount() {
        this.subscribe(this.plugin.state.data.behaviors.isUpdating, v => {
            this.setState({ isDisabled: v });
        });

        this.subscribe(this.plugin.state.events.cell.stateUpdated, e => {
            if (!this.state.isDisabled && MesoscaleState.has(this.plugin) && MesoscaleState.ref(this.plugin) === e.ref) {
                this.forceUpdate();
            }
        });

        this.subscribe(this.plugin.managers.structure.selection.events.changed, e => {
            if (!this.state.isDisabled) {
                this.forceUpdate();
            }
        });
    }

    render() {
        return <div className='msp-scrollable-container'>
            {this.hasModelInfo && <>
                <SectionHeader title='Model' />
                <ModelInfo />
                <Spacer />
            </>}

            <>
                <SectionHeader title='Selection' />
                <SelectionInfo />
                <Spacer />
                <StructureMeasurementsControls initiallyCollapsed={true}/>
            </>
            <MesoQuickStylesControls />
            <Spacer />
            <SectionHeader title='Entities' />
            <EntityControls />
            <Spacer />
            {this.hasFocusInfo && <>
                <SectionHeader title='Focus Info' />
                <FocusInfo />
                <Spacer />
            </>}
        </div>;
    }
}
