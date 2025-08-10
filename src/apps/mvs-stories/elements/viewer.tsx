/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolViewSpec } from '../../../extensions/mvs/behavior';
import { loadMVSData } from '../../../extensions/mvs/components/formats';
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
                this.context.state.isLoading.next(true);
                if (cmd.kind === 'load-mvs') {
                    let loadedData: MVSData | StringLike | Uint8Array | undefined;
                    if (cmd.url) {
                        const data = await this.plugin.runTask(this.plugin.fetch({ url: cmd.url, type: cmd.format === 'mvsx' ? 'binary' : 'string' }));
                        loadedData = await loadMVSData(this.plugin, data, cmd.format ?? 'mvsj', { sourceUrl: cmd.url });
                    } else if (cmd.data) {
                        loadedData = await loadMVSData(this.plugin, cmd.data, cmd.format ?? 'mvsj');
                    }
                    if (StringLike.is(loadedData) || loadedData instanceof Uint8Array) {
                        this.context.state.currentStoryData.next(loadedData as string | Uint8Array);
                    } else if (loadedData) {
                        this.context.state.currentStoryData.next(JSON.stringify(loadedData));
                    }
                }
            } catch (e) {
                console.error(e);
                PluginCommands.Toast.Show(
                    this.plugin,
                    { key: '<mvsload>', title: 'Error', message: e?.message ? `${e?.message}` : `${e}`, timeoutMs: 10000 }
                );
            } finally {
                this.context.state.isLoading.next(false);
            }
        });

        const viewers = this.context.state.viewers.value;
        const next = [...viewers, { name: this.options?.name, model: this }];
        this.context.state.viewers.next(next);
    }

    constructor(private options?: { context?: { name?: string, container?: object }, name?: string }) {
        super();

        this.context = getMVSStoriesContext(options?.context);

        const viewers = this.context.state.viewers.value;
        const index = viewers.findIndex(v => v.name === options?.name);
        if (index >= 0) {
            const next = [...viewers];
            next[index].model.dispose();
            next.splice(index, 0);
            this.context.state.viewers.next(next);
        }
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