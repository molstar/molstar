/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mp4EncoderUI } from '../../../extensions/mp4-export/ui';
import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { AnimationViewportControls, LociLabels, SelectionViewportControls, StateSnapshotViewportControls, TrajectoryViewportControls, ViewportSnapshotDescription } from '../../../mol-plugin-ui/controls';
import { SectionHeader } from '../../../mol-plugin-ui/controls/common';
import { StructureMeasurementsControls } from '../../../mol-plugin-ui/structure/measurements';
import { BackgroundTaskProgress } from '../../../mol-plugin-ui/task';
import { Toasts } from '../../../mol-plugin-ui/toast';
import { Viewport, ViewportControls } from '../../../mol-plugin-ui/viewport';
import { MesoscaleExplorerState } from '../app';
import { MesoscaleState } from '../data/state';
import { CanvasInfo, EntityControls, ModelInfo, SelectionInfo } from './entities';
import { LoaderControls, ExampleControls, SessionControls, SnapshotControls, DatabaseControls } from './states';

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
                <ViewportSnapshotDescription />
            </div>
            <SelectionViewportControls />
            <VPControls />
            <BackgroundTaskProgress />
            <div className='msp-highlight-toast-wrapper'>
                <LociLabels />
                <Toasts />
                <CanvasInfo />
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