/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginUIComponent } from 'mol-plugin/ui/base';
import { CurrentObject } from 'mol-plugin/ui/plugin';
import { AnimationControls } from 'mol-plugin/ui/state/animation';
import { CameraSnapshots } from 'mol-plugin/ui/camera';

export class ControlsWrapper extends PluginUIComponent {
    render() {
        return <div className='msp-scrollable-container msp-right-controls'>
            <CurrentObject />
            <AnimationControls />
            <CameraSnapshots />
        </div>;
    }
}