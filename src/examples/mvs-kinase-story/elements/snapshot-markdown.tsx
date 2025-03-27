/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, distinctUntilChanged, map } from 'rxjs';
import { PluginComponent } from '../../../mol-plugin-state/component';
import { getMolComponentContext, MolComponentContext } from '../context';
import { MolComponentViewerModel } from './viewer';
import Markdown from 'react-markdown';
import { useBehavior } from '../../../mol-plugin-ui/hooks/use-behavior';
import { createRoot } from 'react-dom/client';
import { PluginStateSnapshotManager } from '../../../mol-plugin-state/manager/snapshots';
import { MarkdownAnchor } from '../../../mol-plugin-ui/controls';
import { PluginReactContext } from '../../../mol-plugin-ui/base';

export class MolComponentSnapshotMarkdownModel extends PluginComponent {
    readonly context: MolComponentContext;
    root: HTMLElement | undefined = undefined;

    state = new BehaviorSubject<{
        entry?: PluginStateSnapshotManager.Entry,
        index?: number,
        all: PluginStateSnapshotManager.Entry[],
    }>({ all: [] });

    get viewer() {
        return this.context.behavior.viewers.value?.find(v => this.options?.viewerName === v.name);
    }

    sync() {
        const mng = this.viewer?.model.plugin?.managers.snapshot;
        this.state.next({
            entry: mng?.current,
            index: mng?.current ? mng?.getIndex(mng.current) : undefined,
            all: mng?.state.entries.toArray() ?? [],
        });
    }

    async mount(root: HTMLElement) {
        this.root = root;

        createRoot(root).render(<MolComponentSnapshotMarkdownUI model={this} />);

        let currentViewer: MolComponentViewerModel | undefined = undefined;
        let sub: { unsubscribe: () => void } | undefined = undefined;
        this.subscribe(this.context.behavior.viewers.pipe(
            map(xs => xs.find(v => this.options?.viewerName === v.name)),
            distinctUntilChanged((a, b) => a?.model === b?.model)
        ), viewer => {
            if (currentViewer !== viewer) {
                currentViewer = viewer?.model;
                sub?.unsubscribe();
            }
            if (!viewer) return;
            sub = this.subscribe(viewer.model.plugin?.managers.snapshot.events.changed, () => {
                this.sync();
            });
            this.sync();
        });

        this.sync();
    }

    constructor(private options?: { context?: { name?: string, container?: object }, viewerName?: string }) {
        super();

        this.context = getMolComponentContext(options?.context);
    }
}

export function MolComponentSnapshotMarkdownUI({ model }: { model: MolComponentSnapshotMarkdownModel }) {
    const state = useBehavior(model.state);

    if (state.all.length === 0) {
        return <div>
            <i>No snapshot loaded</i>
        </div>;
    }

    return <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
        <div style={{ display: 'flex', flexDirection: 'row', width: '100%', gap: '8px' }} className='mc-snapshot-markdown-header'>
            <span style={{ lineHeight: '38px', minWidth: 60, maxWidth: 60, flexShrink: 0 }}>{typeof state.index === 'number' ? state.index + 1 : '-'}/{state.all.length}</span>
            <button onClick={() => model.viewer?.model.plugin?.managers.snapshot.applyNext(-1)} style={{ flexGrow: 1, flexShrink: 0 }}>Prev</button>
            <button onClick={() => model.viewer?.model.plugin?.managers.snapshot.applyNext(1)} style={{ flexGrow: 1, flexShrink: 0 }}>Next</button>
        </div>
        <div style={{ flexGrow: 1, overflow: 'hidden', overflowY: 'auto', position: 'relative' }}>
            <div style={{ position: 'absolute', inset: 0 }}>
                <PluginReactContext.Provider value={model.viewer?.model.plugin as any}>
                    <Markdown skipHtml components={{ a: MarkdownAnchor }}>{state.entry?.description ?? 'Description not available'}</Markdown>
                </PluginReactContext.Provider>
            </div>
        </div>
    </div>;
}

export class MolComponentSnapshotMarkdownViewer extends HTMLElement {
    private model: MolComponentSnapshotMarkdownModel | undefined = undefined;

    async connectedCallback() {
        this.model = new MolComponentSnapshotMarkdownModel({
            context: { name: this.getAttribute('context-name') ?? undefined },
            viewerName: this.getAttribute('viewer-name') ?? undefined,
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

window.customElements.define('mc-snapshot-markdown', MolComponentSnapshotMarkdownViewer);