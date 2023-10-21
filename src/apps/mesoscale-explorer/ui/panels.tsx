/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { SectionHeader } from '../../../mol-plugin-ui/controls/common';
import { MesoscaleExplorerState } from '../app';
import { MesoscaleState } from '../data/state';
import { EntityControls, ModelInfo } from './entities';
import { LoaderControls, ExampleControls, SessionControls, SnapshotControls, DatabaseControls } from './states';
import { debounceTime, filter } from 'rxjs/operators';

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

export class RightPanel extends PluginUIComponent {
    get hasInfo() {
        return (
            MesoscaleState.has(this.plugin) &&
            (MesoscaleState.get(this.plugin).description ||
                MesoscaleState.get(this.plugin).link)
        );
    }

    componentDidMount() {
        this.subscribe(this.plugin.state.events.cell.stateUpdated.pipe(filter(e => (MesoscaleState.has(this.plugin) && MesoscaleState.ref(this.plugin) === e.ref)), debounceTime(33)), e => {
            this.forceUpdate();
        });
    }

    render() {
        return <div className='msp-scrollable-container'>
            {this.hasInfo && <>
                <SectionHeader title='Model' />
                <ModelInfo />
                <Spacer />
            </>}

            <SectionHeader title='Entities' />
            <EntityControls />
        </div>;
    }
}