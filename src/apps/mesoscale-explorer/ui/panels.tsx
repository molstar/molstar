/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { SectionHeader } from '../../../mol-plugin-ui/controls/common';
import { MesoscaleExplorerState } from '../app';
import { EntityControls } from './entities';
import { LoaderControls, ExampleControls, SessionControls, SnapshotControls, DatabaseControls } from './states';

export class LeftPanel extends PluginUIComponent {
    render() {
        const customState = this.plugin.customState as MesoscaleExplorerState;

        return <div className='msp-scrollable-container'>
            <SectionHeader title='Database' />
            <DatabaseControls />

            <SectionHeader title='Model' />
            <LoaderControls />

            {customState.examples?.length && <>
                <SectionHeader title='Example' />
                <ExampleControls />
            </>}

            <SectionHeader title='Session' />
            <SessionControls />

            <SectionHeader title='Snapshots' />
            <SnapshotControls />
        </div>;
    }
}

export class RightPanel extends PluginUIComponent {
    render() {
        return <div className='msp-scrollable-container'>
            <SectionHeader title='Entities' />
            <EntityControls />
        </div>;
    }
}