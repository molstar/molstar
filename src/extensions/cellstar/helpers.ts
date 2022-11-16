import { Volume } from '../../mol-model/volume';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { setSubtreeVisibility } from '../../mol-plugin/behavior/static/state';
import { StateBuilder, StateObjectSelector, StateTransformer } from '../../mol-state';
import { ParamDefinition } from '../../mol-util/param-definition';


/** Split entry ID (e.g. 'emd-1832') into source ('emdb') and number ('1832') */
export function splitEntryId(entryId: string) {
    const PREFIX_TO_SOURCE: { [prefix: string]: string } = { 'empiar': 'empiar', 'emd': 'emdb' };
    const [prefix, entry] = entryId.split('-');
    return {
        source: PREFIX_TO_SOURCE[prefix],
        entryNumber: entry
    };
}

/** Create entry ID (e.g. 'emd-1832') for a combination of source ('emdb') and number ('1832') */
export function createEntryId(source: string, entryNumber: string | number) {
    const SOURCE_TO_PREFIX: { [prefix: string]: string } = { 'empiar': 'empiar', 'emdb': 'emd' };
    return `${SOURCE_TO_PREFIX[source]}-${entryNumber}`;
}



/**
 * Represents a set of values to choose from, with a default value. Example:
 * ```
 * export const MyChoice = new Choice({ yes: 'I agree', no: 'Nope' }, 'yes');
 * export type MyChoiceType = Choice.Values<typeof MyChoice>; // 'yes'|'no'
 * ```
 */
export class Choice<T extends string, D extends T> {
    readonly defaultValue: D;
    readonly options: [T, string][];
    private readonly nameDict: { [value in T]: string };
    constructor(opts: { [value in T]: string }, defaultValue: D) {
        this.defaultValue = defaultValue;
        this.options = Object.keys(opts).map(k => [k as T, opts[k as T]]);
        this.nameDict = opts;
    }
    PDSelect(defaultValue?: T, info?: ParamDefinition.Info): ParamDefinition.Select<T> {
        return ParamDefinition.Select<T>(defaultValue ?? this.defaultValue, this.options, info);
    }
    prettyName(value: T): string {
        return this.nameDict[value];
    }
}
export namespace Choice {
    export type Values<T extends Choice<any, any>> = T extends Choice<infer R, any> ? R : any;
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
            setSubtreeVisibility(node.state!, node.ref, true);  // hide
        }
    }

    public async showNode(key: string, factory: () => StateObjectSelector | Promise<StateObjectSelector>, forceVisible: boolean = true) {
        console.log('showNode:', key);
        let node = this.getNode(key);
        if (node) {
            console.log('showNode set visible');
            if (forceVisible) {
                setSubtreeVisibility(node.state!, node.ref, false);  // show
            }
        } else {
            console.log('showNode create');
            node = await factory();
            this.nodes[key] = node;
        }
        return node;
    }

}


const CreateTransformer = StateTransformer.builderFactory('cellstar');

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
})
