/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolViewSpec } from '../../../extensions/mvs/behavior';
import { loadMVSX } from '../../../extensions/mvs/components/formats';
import { loadMVS } from '../../../extensions/mvs/load';
import { MVSData } from '../../../extensions/mvs/mvs-data';
import { StringLike } from '../../../mol-io/common/string-like';
import { PluginComponent } from '../../../mol-plugin-state/component';
import { createPluginUI } from '../../../mol-plugin-ui';
import { renderReact18 } from '../../../mol-plugin-ui/react18';
import { DefaultPluginUISpec } from '../../../mol-plugin-ui/spec';
import { PluginCommands } from '../../../mol-plugin/commands';
import { PluginConfig } from '../../../mol-plugin/config';
import { PluginContext } from '../../../mol-plugin/context';
import { PluginSpec } from '../../../mol-plugin/spec';
import { Task } from '../../../mol-task';
import { getMVSStoriesContext, MVSStoriesContext } from '../context';

export class MVSStoriesViewerModel extends PluginComponent {
    readonly context: MVSStoriesContext;
    plugin?: PluginContext = undefined;

    async mount(root: HTMLElement) {
        const spec = DefaultPluginUISpec();
        this.plugin = await createPluginUI({
            target: root,
            render: renderReact18,
            spec: {
                ...spec,
                layout: {
                    initial: {
                        isExpanded: false,
                        showControls: false,
                        controlsDisplay: 'landscape',
                    },
                },
                components: {
                    remoteState: 'none',
                    viewport: {
                        snapshotDescription: EmptyDescription,
                    }
                },
                behaviors: [
                    ...spec.behaviors,
                    PluginSpec.Behavior(MolViewSpec)
                ],
                config: [
                    [PluginConfig.Viewport.ShowAnimation, false],
                ]
            }
        });

        this.subscribe(this.context.commands, async (cmd) => {
            if (!cmd || !this.plugin) return;

            try {
                if (cmd.kind === 'load-mvs') {
                    if (cmd.url) {
                        const data = await this.plugin.runTask(this.plugin!.fetch({ url: cmd.url, type: cmd.format === 'mvsx' ? 'binary' : 'string' }));
                        await loadMvsData(this.plugin, data, cmd.format ?? 'mvsj', { sourceUrl: cmd.url });
                    } else if (cmd.data) {
                        await loadMvsData(this.plugin, cmd.data, cmd.format ?? 'mvsj');
                    }
                }
            } catch (e) {
                console.error(e);
                PluginCommands.Toast.Show(
                    this.plugin,
                    { key: '<mvsload>', title: 'Error', message: e?.message ? `${e?.message}` : `${e}`, timeoutMs: 10000 }
                );
            }
        });

        const viewers = this.context.behavior.viewers.value;
        const next = [...viewers, { name: this.options?.name, model: this }];
        this.context.behavior.viewers.next(next);
    }

    constructor(private options?: { context?: { name?: string, container?: object }, name?: string }) {
        super();

        this.context = getMVSStoriesContext(options?.context);

        const viewers = this.context.behavior.viewers.value;
        const index = viewers.findIndex(v => v.name === options?.name);
        if (index >= 0) {
            const next = [...viewers];
            next[index].model.dispose();
            next.splice(index, 0);
            this.context.behavior.viewers.next(next);
        }
    }
}

async function loadMvsData(plugin: PluginContext, data: MVSData | StringLike | Uint8Array, format: 'mvsj' | 'mvsx', options?: { sourceUrl?: string }) {
    if (typeof data === 'string' && data.startsWith('base64')) {
        data = Uint8Array.from(atob(data.substring(7)), c => c.charCodeAt(0)); // Decode base64 string to Uint8Array
    }

    if (format === 'mvsj') {
        if ((data as Uint8Array).BYTES_PER_ELEMENT && (data as Uint8Array).buffer) {
            data = new TextDecoder().decode(data as Uint8Array); // Decode Uint8Array to string using UTF8
        }

        let mvsData: MVSData;
        if (typeof data === 'string') {
            mvsData = MVSData.fromMVSJ(data);
        } else {
            mvsData = data as MVSData;
        }
        await loadMVS(plugin, mvsData, { sanityChecks: true, sourceUrl: undefined, ...options });
    } else if (format === 'mvsx') {
        if (typeof data === 'string') {
            throw new Error("loadMvsData: if `format` is 'mvsx', then `data` must be a Uint8Array or a base64-encoded string prefixed with 'base64,'.");
        }
        await plugin.runTask(Task.create('Load MVSX file', async ctx => {
            const parsed = await loadMVSX(plugin, ctx, data as Uint8Array);
            await loadMVS(plugin, parsed.mvsData, { sanityChecks: true, sourceUrl: parsed.sourceUrl, ...options });
        }));
    } else {
        throw new Error(`Unknown MolViewSpec format: ${format}`);
    }
}

function EmptyDescription() {
    return <></>;
}

export class MVSStoriesViewer extends HTMLElement {
    private model: MVSStoriesViewerModel | undefined = undefined;

    async connectedCallback() {
        this.model = new MVSStoriesViewerModel({
            name: this.getAttribute('name') ?? undefined,
            context: { name: this.getAttribute('context-name') ?? undefined },
        });
        await this.model.mount(this);
    }

    disconnectedCallback() {
        this.model?.dispose();
        this.model = undefined;
    }

    constructor() {
        super();
    }
}

window.customElements.define('mvs-stories-viewer', MVSStoriesViewer);