/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as ReactDOM from 'react-dom';
import { PluginUIContext } from '@molstar/plugin-ui/context';
import { PluginContextContainer } from '@molstar/plugin-ui/plugin';
import { TransformUpdaterControl } from '@molstar/plugin-ui/state/update-transform';
import { StateElements } from '@molstar/proteopedia-wrapper-example/helpers';

export function volumeStreamingControls(plugin: PluginUIContext, parent: Element) {
  ReactDOM.render(
    <PluginContextContainer plugin={plugin}>
      <TransformUpdaterControl nodeRef={StateElements.VolumeStreaming} />
    </PluginContextContainer>,
    parent,
  );
}
