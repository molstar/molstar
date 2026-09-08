/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Volume Tools — landing page (`index.html`) linking to one page per tool. Each tool runs in
 * its own document, so the two never share a plugin instance or a scene.
 *
 * Start:  npm run dev -- -e volume-tools
 * Serve:  http-server -p 1338 -g
 * Open:   http://localhost:1338/build/examples/volume-tools/
 *         (a tool page also takes ?url=<map.mrc> to load a volume right away)
 */

import * as React from 'react';
import * as ReactDOM from 'react-dom';
import { createPluginUI } from '../../mol-plugin-ui';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { PluginSpec } from '../../mol-plugin/spec';
import { VolumeMaskBehavior } from '../../extensions/volume-mask';
import { VolumeSegmentorBehavior, VolumeSegmentorManager } from '../../extensions/volume-segmentor';
import { VolumeMaskController } from './controller';
import { BodiesPanel } from './ui/bodies-panel';
import { MaskCreatorPanel } from './ui/mask-panel';
import '../../mol-plugin-ui/skin/light.scss';
import './index.html';
import './tool.html';

interface Tool {
    name: string;
    behavior: PluginSpec.Behavior;
    /** Builds the panel once the plugin is up, and returns anything worth exposing on `window`. */
    panel: (plugin: PluginContext) => { element: React.ReactElement, globals: Record<string, unknown> };
}

const Tools: { [id: string]: Tool } = {
    mask: {
        name: 'Mask Creator',
        behavior: PluginSpec.Behavior(VolumeMaskBehavior),
        panel: plugin => {
            const controller = new VolumeMaskController(plugin);
            return {
                element: React.createElement(MaskCreatorPanel, { plugin, controller }),
                globals: { maskController: controller },
            };
        },
    },
    segmentor: {
        name: 'Segmentor',
        behavior: PluginSpec.Behavior(VolumeSegmentorBehavior),
        panel: plugin => {
            const manager = VolumeSegmentorManager.get(plugin)!;
            return {
                element: React.createElement(BodiesPanel, { plugin, manager }),
                globals: { segmentorManager: manager },
            };
        },
    },
};

async function init() {
    const params = new URL(window.location.href).searchParams;
    const tool = Tools[params.get('tool') ?? ''];
    if (!tool) {
        window.location.replace('./');
        return;
    }

    document.title = `${tool.name} — Volume Tools`;
    document.getElementById('tool-name')!.textContent = tool.name;

    const spec = DefaultPluginUISpec();
    const plugin = await createPluginUI({
        target: document.getElementById('app')!,
        render: renderReact18,
        spec: {
            ...spec,
            behaviors: [...spec.behaviors, tool.behavior],
            layout: {
                initial: {
                    isExpanded: false, // stay inside #app; expanded mode is position: fixed
                    showControls: false, // hide Mol* built-in panels; we use our own
                },
            },
            components: {
                ...spec.components,
                viewport: { ...spec.components?.viewport, controls: () => null },
            },
            config: [
                ...(spec.config ?? []),
                [PluginConfig.Viewport.ShowAnimation, false],
                [PluginConfig.Viewport.ShowTrajectoryControls, false],
            ],
        },
    });

    const { element, globals } = tool.panel(plugin);
    ReactDOM.render(element, document.getElementById('panel-root')!);

    (window as any).plugin = plugin;
    for (const [key, value] of Object.entries(globals)) (window as any)[key] = value;

    const url = params.get('url');
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
