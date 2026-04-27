/**
 * Volume Mask Creator — standalone example page.
 *
 * Start:  npm run dev -- -e volume-mask
 * Serve:  http-server -p 1338
 * Open:   http://localhost:1338/build/examples/volume-mask/
 */

import * as React from 'react';
import * as ReactDOM from 'react-dom';
import { createPluginUI } from '../../mol-plugin-ui';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginSpec } from '../../mol-plugin/spec';
import { VolumeMaskBehavior, VolumeMaskController, MaskCreatorPanel } from '../../extensions/volume-mask';
import '../../mol-plugin-ui/skin/light.scss';
import './index.html';

async function init() {
    const plugin = await createPluginUI({
        target: document.getElementById('app')!,
        render: renderReact18,
        spec: {
            ...DefaultPluginUISpec(),
            behaviors: [
                ...DefaultPluginUISpec().behaviors,
                PluginSpec.Behavior(VolumeMaskBehavior),
            ],
            layout: {
                initial: {
                    isExpanded: true,
                    showControls: false, // hide Mol* built-in panels; we use our own
                },
            },
        },
    });

    const controller = new VolumeMaskController(plugin);
    ReactDOM.render(
        React.createElement(MaskCreatorPanel, { plugin, controller }),
        document.getElementById('mask-panel')!
    );

    (window as any).plugin = plugin;
    (window as any).maskController = controller;
    console.log('Ready. Use the "Open Volume" button to load a local .mrc/.map file.');
}

init().catch(console.error);
