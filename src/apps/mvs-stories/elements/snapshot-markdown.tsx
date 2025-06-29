/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, distinctUntilChanged, map } from 'rxjs';
import { PluginComponent } from '../../../mol-plugin-state/component';
import { getMVSStoriesContext, MVSStoriesContext } from '../context';
import { MVSStoriesViewerModel } from './viewer';
import { useBehavior } from '../../../mol-plugin-ui/hooks/use-behavior';
import { createRoot } from 'react-dom/client';
import { PluginStateSnapshotManager } from '../../../mol-plugin-state/manager/snapshots';
import { PluginReactContext } from '../../../mol-plugin-ui/base';
import { CSSProperties } from 'react';
import { Markdown } from '../../../mol-plugin-ui/controls/markdown';

export class MVSStoriesSnapshotMarkdownModel extends PluginComponent {
    readonly context: MVSStoriesContext;
    root: HTMLElement | undefined = undefined;

    state = new BehaviorSubject<{
        entry?: PluginStateSnapshotManager.Entry,
        index?: number,
        all: PluginStateSnapshotManager.Entry[],
    }>({ all: [] });

    get viewer() {
        return this.context.state.viewers.value?.find(v => this.options?.viewerName === v.name);
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

        createRoot(root).render(<MVSStoriesSnapshotMarkdownUI model={this} />);

        let currentViewer: MVSStoriesViewerModel | undefined = undefined;
        let sub: { unsubscribe: () => void } | undefined = undefined;
        this.subscribe(this.context.state.viewers.pipe(
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

        this.context = getMVSStoriesContext(options?.context);
    }
}

export function MVSStoriesSnapshotMarkdownUI({ model }: { model: MVSStoriesSnapshotMarkdownModel }) {
    const state = useBehavior(model.state);
    const isLoading = useBehavior(model.context.state.isLoading);

    const style: CSSProperties = { display: 'flex', flexDirection: 'column', height: '100%' };
    const className = 'mvs-stories-markdown-explanation';

    if (isLoading) {
        return <div style={style} className={className}>
            <i>Loading...</i>
        </div>;
    }

    if (state.all.length === 0) {
        return <div style={style} className={className}>
            <i>No snapshot loaded or no description available</i>
        </div>;
    }

    return <div style={style} className={className}>
        <div style={{ display: 'flex', flexDirection: 'row', width: '100%', gap: '8px' }}>
            <span style={{ lineHeight: '38px', minWidth: 60, maxWidth: 60, flexShrink: 0 }}>{typeof state.index === 'number' ? state.index + 1 : '-'}/{state.all.length}</span>
            <button onClick={() => model.viewer?.model.plugin?.managers.snapshot.applyNext(-1)} style={{ flexGrow: 1, flexShrink: 0 }}>Prev</button>
            <button onClick={() => model.viewer?.model.plugin?.managers.snapshot.applyNext(1)} style={{ flexGrow: 1, flexShrink: 0 }}>Next</button>
        </div>
        <div style={{ flexGrow: 1, overflow: 'hidden', overflowY: 'auto', position: 'relative' }}>
            <div style={{ position: 'absolute', inset: 0 }}>
                <PluginReactContext.Provider value={model.viewer?.model.plugin as any}>
                    <Markdown>{state.entry?.description ?? 'Description not available'}</Markdown>
                </PluginReactContext.Provider>
            </div>
        </div>
    </div>;
}

export class MVSStoriesSnapshotMarkdownViewer extends HTMLElement {
    private model: MVSStoriesSnapshotMarkdownModel | undefined = undefined;

    async connectedCallback() {
        this.model = new MVSStoriesSnapshotMarkdownModel({
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

window.customElements.define('mvs-stories-snapshot-markdown', MVSStoriesSnapshotMarkdownViewer);