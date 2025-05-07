/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolViewSpec } from '../../../extensions/mvs/behavior';
import { loadMVS } from '../../../extensions/mvs/load';
import { MVSData } from '../../../extensions/mvs/mvs-data';
import { StringLike } from '../../../mol-io/common/string-like';
import { PluginComponent } from '../../../mol-plugin-state/component';
import { createPluginUI } from '../../../mol-plugin-ui';
import { renderReact18 } from '../../../mol-plugin-ui/react18';
import { DefaultPluginUISpec } from '../../../mol-plugin-ui/spec';
import { PluginConfig } from '../../../mol-plugin/config';
import { PluginContext } from '../../../mol-plugin/context';
import { PluginSpec } from '../../../mol-plugin/spec';
import { getMolComponentContext, MolComponentContext } from '../context';

export class MolComponentViewerModel extends PluginComponent {
    readonly context: MolComponentContext;
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
            if (!cmd) return;

            if (cmd.kind === 'load-mvs') {
                if (cmd.url) {
                    const data = await this.plugin!.runTask(this.plugin!.fetch({ url: cmd.url, type: 'string' }));
                    const mvsData = MVSData.fromMVSJ(StringLike.toString(data));
                    await loadMVS(this.plugin!, mvsData, { sanityChecks: true, sourceUrl: cmd.url, replaceExisting: true });
                } else if (cmd.data) {
                    await loadMVS(this.plugin!, cmd.data, { sanityChecks: true, replaceExisting: true });
                }
            }
        });

        const viewers = this.context.behavior.viewers.value;
        const next = [...viewers, { name: this.options?.name, model: this }];
        this.context.behavior.viewers.next(next);
    }

    constructor(private options?: { context?: { name?: string, container?: object }, name?: string }) {
        super();

        this.context = getMolComponentContext(options?.context);

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

function EmptyDescription() {
    return <></>;
}

export class MolComponentViewer extends HTMLElement {
    private model: MolComponentViewerModel | undefined = undefined;

    async connectedCallback() {
        this.model = new MolComponentViewerModel({
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

window.customElements.define('mc-viewer', MolComponentViewer);