/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Volume } from '../../mol-model/volume';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { setSubtreeVisibility } from '../../mol-plugin/behavior/static/state';
import { StateBuilder, StateObjectSelector, StateTransformer } from '../../mol-state';
import { ParamDefinition } from '../../mol-util/param-definition';
import { Source } from './entry-root';


/** Split entry ID (e.g. 'emd-1832') into source ('emdb') and number ('1832') */
export function splitEntryId(entryId: string) {
    const PREFIX_TO_SOURCE: { [prefix: string]: Source } = { 'emd': 'emdb' };
    const [prefix, entry] = entryId.split('-');
    return {
        source: PREFIX_TO_SOURCE[prefix] ?? prefix,
        entryNumber: entry
    };
}

/** Create entry ID (e.g. 'emd-1832') for a combination of source ('emdb') and number ('1832') */
export function createEntryId(source: Source, entryNumber: string | number) {
    const SOURCE_TO_PREFIX: { [prefix: string]: string } = { 'emdb': 'emd' };
    const prefix = SOURCE_TO_PREFIX[source] ?? source;
    return `${prefix}-${entryNumber}`;
}


export function isDefined<T>(x: T | undefined): x is T {
    return x !== undefined;
}


export class NodeManager {
    private nodes: { [key: string]: StateObjectSelector };

    constructor() {
        this.nodes = {};
    }

    private static nodeExists(node: StateObjectSelector): boolean {
        try {
            return node.checkValid();
        } catch {
            return false;
        }
    }

    public getNode(key: string): StateObjectSelector | undefined {
        const node = this.nodes[key];
        if (node && !NodeManager.nodeExists(node)) {
            delete this.nodes[key];
            return undefined;
        }
        return node;
    }

    public getNodes(): StateObjectSelector[] {
        return Object.keys(this.nodes).map(key => this.getNode(key)).filter(node => node) as StateObjectSelector[];
    }

    public deleteAllNodes(update: StateBuilder.Root) {
        for (const node of this.getNodes()) {
            update.delete(node);
        }
        this.nodes = {};
    }

    public hideAllNodes() {
        for (const node of this.getNodes()) {
            setSubtreeVisibility(node.state!, node.ref, true); // hide
        }
    }

    public async showNode(key: string, factory: () => StateObjectSelector | Promise<StateObjectSelector>, forceVisible: boolean = true) {
        let node = this.getNode(key);
        if (node) {
            if (forceVisible) {
                setSubtreeVisibility(node.state!, node.ref, false); // show
            }
        } else {
            node = await factory();
            this.nodes[key] = node;
        }
        return node;
    }
}



const CreateTransformer = StateTransformer.builderFactory('volseg');

export const CreateVolume = CreateTransformer({
    name: 'create-transformer',
    from: PluginStateObject.Root,
    to: PluginStateObject.Volume.Data,
    params: {
        label: ParamDefinition.Text('Volume', { isHidden: true }),
        description: ParamDefinition.Text('', { isHidden: true }),
        volume: ParamDefinition.Value<Volume>(undefined as any, { isHidden: true }),
    }
})({
    apply({ params }) {
        return new PluginStateObject.Volume.Data(params.volume, { label: params.label, description: params.description });
    }
});



export function applyEllipsis(name: string, max_chars: number = 60) {
    if (name.length <= max_chars) return name;
    const beginning = name.substring(0, max_chars);
    let lastSpace = beginning.lastIndexOf(' ');
    if (lastSpace === -1) return beginning + '...';
    if (lastSpace > 0 && ',;.'.includes(name.charAt(lastSpace - 1))) lastSpace--;
    return name.substring(0, lastSpace) + '...';
}


export function lazyGetter<T>(getter: () => T, errorIfUndefined?: string) {
    let value: T | undefined = undefined;
    return () => {
        if (value === undefined) value = getter();
        if (errorIfUndefined && value === undefined) throw new Error(errorIfUndefined);
        return value;
    };
}
