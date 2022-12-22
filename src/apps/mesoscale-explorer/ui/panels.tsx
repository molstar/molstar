/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { EntityControls } from './entities';
import { LoaderControls } from './states';

export class LeftPanel extends PluginUIComponent {
    render() {
        return <div className='msp-scrollable-container'>
            <LoaderControls />
        </div>;
    }
}

export class RightPanel extends PluginUIComponent {
    render() {
        return <div className='msp-scrollable-container'>
            <EntityControls />
        </div>;
    }
}