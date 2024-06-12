/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { throttleTime } from 'rxjs';
import { Mp4EncoderUI } from '../../../extensions/mp4-export/ui';
import { CollapsableControls, CollapsableState, PluginUIComponent } from '../../../mol-plugin-ui/base';
import { AnimationViewportControls, LociLabels, SelectionViewportControls, StateSnapshotViewportControls, TrajectoryViewportControls } from '../../../mol-plugin-ui/controls';
import { SectionHeader } from '../../../mol-plugin-ui/controls/common';
import { TuneSvg } from '../../../mol-plugin-ui/controls/icons';
import { StructureMeasurementsControls } from '../../../mol-plugin-ui/structure/measurements';
import { BackgroundTaskProgress } from '../../../mol-plugin-ui/task';
import { Toasts } from '../../../mol-plugin-ui/toast';
import { Viewport, ViewportControls } from '../../../mol-plugin-ui/viewport';
import { PluginCommands } from '../../../mol-plugin/commands';
import { MesoscaleExplorerState } from '../app';
import { MesoscaleState } from '../data/state';
import { EntityControls, ModelInfo, SelectionInfo } from './entities';
import { LoaderControls, ExampleControls, SessionControls, SnapshotControls, DatabaseControls, MesoViewportSnapshotDescription } from './states';
import { ParameterControls } from '../../../mol-plugin-ui/controls/parameters';
import { Canvas3DContext, Canvas3DParams } from '../../../mol-canvas3d/canvas3d';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

const Spacer = () => <div style={{ height: '2em' }} />;

export class MesoScaleViewport extends PluginUIComponent {
    render() {
        const VPControls = this.plugin.spec.components?.viewport?.controls || ViewportControls;
        return <>
            <Viewport />
            <div className='msp-viewport-top-left-controls'>
                <AnimationViewportControls />
                <TrajectoryViewportControls />
                <StateSnapshotViewportControls />
                <MesoViewportSnapshotDescription />
            </div>
            <SelectionViewportControls />
            <VPControls />
            <BackgroundTaskProgress />
            <div className='msp-highlight-toast-wrapper'>
                <LociLabels />
                <Toasts />
            </div>
        </>;
    }
}

export class LeftPanel extends PluginUIComponent {
    render() {
        const customState = this.plugin.customState as MesoscaleExplorerState;

        return <div className='msp-scrollable-container'>
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
                <StructureMeasurementsControls />
            </>

            <SectionHeader title='Entities' />
            <EntityControls />
        </div>;
    }
}


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
                <SectionHeader title='Viewport' />
                <ParameterControls params={Canvas3DParams} values={this.plugin.canvas3d.props} onChange={this.setSettings} />
                <ParameterControls params={Canvas3DContext.Params} values={this.plugin.canvas3dContext.props} onChange={this.setCanvas3DContextProps} />
            </>}
        </>; // Add closing tag for the JSX element
    }

    private setSettings = (p: { param: PD.Base<any>, name: string, value: any }) => {
        PluginCommands.Canvas3D.SetSettings(this.plugin, { settings: { [p.name]: p.value } });
    };

    private setCanvas3DContextProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        this.plugin.canvas3dContext?.setProps({ [p.name]: p.value });
        this.plugin.events.canvas3d.settingsUpdated.next(void 0);
    };

    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
        this.subscribe(this.plugin.layout.events.updated, () => this.forceUpdate());

        if (this.plugin.canvas3d) {
            this.subscribe(this.plugin.canvas3d.camera.stateChanged.pipe(throttleTime(500, undefined, { leading: true, trailing: true })), state => {
                if (state.radiusMax !== undefined || state.radius !== undefined) {
                    this.forceUpdate();
                }
            });
        }
    }
}
