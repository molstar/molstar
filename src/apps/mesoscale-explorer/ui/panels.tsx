/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../../../mol-model/loci';
import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { SectionHeader } from '../../../mol-plugin-ui/controls/common';
import { MesoscaleExplorerState } from '../app';
import { MesoscaleState } from '../data/state';
import { EntityControls, ModelInfo, SelectionInfo } from './entities';
import { LoaderControls, ExampleControls, SessionControls, SnapshotControls, DatabaseControls } from './states';

const Spacer = () => <div style={{ height: '2em' }} />;

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

    get hasSelectionInfo() {
        let count = 0;
        this.plugin.managers.structure.selection.entries.forEach(e => {
            if (!Loci.isEmpty(e.selection)) count += 1;
        });
        return count > 0;
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

            {this.hasSelectionInfo && <>
                <SectionHeader title='Selection' />
                <SelectionInfo />
                <Spacer />
            </>}

            <SectionHeader title='Entities' />
            <EntityControls />
        </div>;
    }
}