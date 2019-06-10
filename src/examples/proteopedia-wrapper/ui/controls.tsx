/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import * as ReactDOM from 'react-dom';
import { PluginUIComponent } from '../../../mol-plugin/ui/base';
import { CurrentObject, PluginContextContainer } from '../../../mol-plugin/ui/plugin';
import { AnimationControls } from '../../../mol-plugin/ui/state/animation';
import { CameraSnapshots } from '../../../mol-plugin/ui/camera';
import { PluginContext } from '../../../mol-plugin/context';
import { TransformUpdaterControl } from '../../../mol-plugin/ui/state/update-transform';
import { StateElements } from '../helpers';

export class ControlsWrapper extends PluginUIComponent {
    render() {
        return <div className='msp-scrollable-container msp-right-controls'>
            <CurrentObject />
            <AnimationControls />
            <CameraSnapshots />
        </div>;
    }
}

export function volumeStreamingControls(plugin: PluginContext, parent: Element) {
    ReactDOM.render(<PluginContextContainer plugin={plugin}>
            <TransformUpdaterControl nodeRef={StateElements.VolumeStreaming} />
        </PluginContextContainer>,
        parent);
}