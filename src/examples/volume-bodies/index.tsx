/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Volume Bodies — standalone example page.
 *
 * Start:  npm run dev -- -e volume-bodies
 * Serve:  http-server -p 1338 -g
 * Open:   http://localhost:1338/build/examples/volume-bodies/ (optionally ?url=<map.mrc>)
 */

import * as React from 'react';
import * as ReactDOM from 'react-dom';
import { createPluginUI } from '../../mol-plugin-ui';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginContext } from '../../mol-plugin/context';
import { PluginSpec } from '../../mol-plugin/spec';
import { VolumeBodiesBehavior, VolumeBodiesManager } from '../../extensions/volume-bodies';
import { BodiesPanel } from './ui/bodies-panel';
import '../../mol-plugin-ui/skin/light.scss';
import './index.html';

async function init() {
    const spec = DefaultPluginUISpec();
    const plugin = await createPluginUI({
        target: document.getElementById('app')!,
        render: renderReact18,
        spec: {
            ...spec,
            behaviors: [...spec.behaviors, PluginSpec.Behavior(VolumeBodiesBehavior)],
            layout: {
                initial: {
                    isExpanded: false, // stay inside #app; expanded mode is position: fixed
                    showControls: false, // hide Mol* built-in panels; we use our own
                },
            },
        },
    });

    const manager = VolumeBodiesManager.get(plugin)!;
    ReactDOM.render(
        React.createElement(BodiesPanel, { plugin, manager }),
        document.getElementById('bodies-panel')!
    );

    (window as any).plugin = plugin;
    (window as any).bodiesManager = manager;

    // Optional: ?url=<volume.mrc> loads a volume right away.
    const url = new URL(window.location.href).searchParams.get('url');
    if (url) await loadVolume(plugin, url);
    console.log('Ready. Use the "Open Volume" button to load a local .mrc/.map file.');
}

async function loadVolume(plugin: PluginContext, url: string) {
    const data = await plugin.builders.data.download({ url, isBinary: true, label: url.split('/').pop() }, { state: { isGhost: true } });
    const provider = plugin.dataFormats.get('ccp4');
    if (!provider) return;
    const parsed = await provider.parse(plugin, data);
    await provider.visuals?.(plugin, parsed);
}

init().catch(console.error);
