/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as ReactDOM from 'react-dom';
import { PluginUIContext } from '../../../mol-plugin-ui/context';
import { PluginContextContainer } from '../../../mol-plugin-ui/plugin';
import { TransformUpdaterControl } from '../../../mol-plugin-ui/state/update-transform';
import { StateElements } from '../helpers';

export function volumeStreamingControls(plugin: PluginUIContext, parent: Element) {
    ReactDOM.render(<PluginContextContainer plugin={plugin}>
        <TransformUpdaterControl nodeRef={StateElements.VolumeStreaming} />
    </PluginContextContainer>, parent);
}