/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { SectionHeader } from '../../../mol-plugin-ui/controls/common';
import { EntityControls } from './entities';
import { LoaderControls, SessionControls, SnapshotControls } from './states';

export class LeftPanel extends PluginUIComponent {
    render() {
        return <div className='msp-scrollable-container'>
            <SectionHeader title='Model' />
            <LoaderControls />

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