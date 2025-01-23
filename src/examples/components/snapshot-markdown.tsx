/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, distinctUntilChanged, map } from 'rxjs';
import { PluginComponent } from '../../mol-plugin-state/component';
import { getMolComponentContext, MolComponentContext } from './context';
import { MolComponentViewerModel } from './viewer';
import Markdown from 'react-markdown';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';
import { createRoot } from 'react-dom/client';

export class MolComponentSnapshotMarkdownModel extends PluginComponent {
    readonly context: MolComponentContext;
    root: HTMLElement | undefined = undefined;

    state = new BehaviorSubject<string>('');

    get viewer() {
        return this.context.behavior.viewers.value?.find(v => this.options?.viewerName === v.name);
    }

    sync() {
        const snapshot = this.viewer?.model.plugin?.managers.snapshot.current;
        this.state.next(snapshot?.description ?? 'no snapshot');
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
    return <div>
        <button onClick={() => model.viewer?.model.plugin?.managers.snapshot.applyNext(-1)}>Prev</button>
        <button onClick={() => model.viewer?.model.plugin?.managers.snapshot.applyNext(1)}>Next</button>
        <Markdown>{state}</Markdown>
    </div>;
}

export class MolComponentSnapshotMarkdownViewer extends HTMLElement {
    async connectedCallback() {
        const model = new MolComponentSnapshotMarkdownModel();
        await model.mount(this);
    }
    constructor() {
        super();
    }
}

window.customElements.define('mol-component-snapshot-markdown', MolComponentSnapshotMarkdownViewer);