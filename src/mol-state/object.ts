/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from '../mol-util';
import { StateTransform } from './transform';
import { ParamDefinition } from '../mol-util/param-definition';
import { State } from './state';
import { StateSelection, StateTransformer } from '../mol-state';
import { StateBuilder } from './state/builder';

export { StateObject, StateObjectCell };

interface StateObject<D = any, T extends StateObject.Type = StateObject.Type<any>> {
    readonly id: UUID,
    readonly type: T,
    readonly data: D,
    readonly label: string,
    readonly description?: string,
    // assigned by reconciler to be StateTransform.props.tag
    readonly tags?: string[]
}

namespace StateObject {
    export function factory<T extends Type>() {
        return <D = { }>(type: T) => create<D, T>(type);
    }

    export type Type<Cls extends string = string> = { name: string, typeClass: Cls }
    export type Ctor<T extends StateObject = StateObject> = { new(...args: any[]): T, is(obj?: StateObject): boolean, type: any }
    export type From<C extends Ctor> = C extends Ctor<infer T> ? T : never

    export function create<Data, T extends Type>(type: T) {
        return class O implements StateObject<Data, T> {
            static type = type;
            static is(obj?: StateObject): obj is O { return !!obj && type === obj.type; }
            id = UUID.create22();
            type = type;
            label: string;
            description?: string;
            constructor(public data: Data, props?: { label: string, description?: string }) {
                this.label = props && props.label || type.name;
                this.description = props && props.description;
            }
        };
    }

    export function hasTag(o: StateObject, t: string) {
        if (!o.tags) return false;
        for (const s of o.tags) {
            if (s === t) return true;
        }
        return false;
    }

    /** A special object indicating a transformer result has no value. */
    export const Null: StateObject<any, any> = {
        id: UUID.create22(),
        type: { name: 'Null', typeClass: 'Null' },
        data: void 0,
        label: 'Null'
    };
}

interface StateObjectCell<T extends StateObject = StateObject, F extends StateTransform = StateTransform> {
    parent?: State,

    transform: F,

    // Which object was used as a parent to create data in this cell
    sourceRef: StateTransform.Ref | undefined,

    status: StateObjectCell.Status,
    state: StateTransform.State,

    params: {
        definition: ParamDefinition.Params,
        values: any
    } | undefined,

    dependencies: {
        dependentBy: StateObjectCell[],
        dependsOn: StateObjectCell[]
    },

    errorText?: string,
    obj?: T,

    cache: unknown | undefined
}

namespace StateObjectCell {
    export type Status = 'ok' | 'error' | 'pending' | 'processing'

    export type Obj<C extends StateObjectCell> = C extends StateObjectCell<infer T> ? T : never
    export type Transform<C extends StateObjectCell> = C extends StateObjectCell<any, infer T> ? T : never
    export type Transformer<C extends StateObjectCell> = C extends StateObjectCell<any, StateTransform<infer T>> ? T : never

    export function is(o: any): o is StateObjectCell {
        const c: StateObjectCell = o;
        return !!c && !!c.transform && !!c.parent && !!c.status;
    }

    export type Ref = StateTransform.Ref | StateObjectCell | StateObjectSelector

    export function resolve(state: State, refOrCellOrSelector: StateTransform.Ref | StateObjectCell | StateObjectSelector) {
        const ref = typeof refOrCellOrSelector === 'string'
            ? refOrCellOrSelector
            : StateObjectCell.is(refOrCellOrSelector)
                ? refOrCellOrSelector.transform.ref
                : refOrCellOrSelector.ref;
        return state.cells.get(ref);
    }
}

// TODO: improve the API?
export class StateObjectTracker<T extends StateObject> {
    private query: StateSelection.Query;
    private version: string = '';
    cell: StateObjectCell | undefined;
    data: T['data'] | undefined;

    setQuery(sel: StateSelection.Selector) {
        this.query = StateSelection.compile(sel);
    }

    update() {
        const cell = this.state.select(this.query)[0];
        const version = cell ? cell.transform.version : void 0;
        const changed = this.cell !== cell || this.version !== version;
        this.cell = cell;
        this.version = version || '';
        this.data = cell && cell.obj ? cell.obj.data as T : void 0 as any;
        return changed;
    }

    constructor(private state: State) { }
}

export class StateObjectSelector<S extends StateObject = StateObject, T extends StateTransformer = StateTransformer> {
    get cell(): StateObjectCell<S, StateTransform<T>> | undefined {
        return this.state?.cells.get(this.ref) as StateObjectCell<S, StateTransform<T>> | undefined;
    }

    get obj(): S | undefined {
        return this.state?.cells.get(this.ref)?.obj as S | undefined;
    }

    get data(): S['data'] | undefined {
        return this.obj?.data;
    }

    /** Create a new build and apply update or use the provided one. */
    update(params: StateTransformer.Params<T>, builder?: StateBuilder.Root | StateBuilder.To<any>): StateBuilder
    update(params: (old: StateTransformer.Params<T>) => StateTransformer.Params<T> | void, builder?: StateBuilder.Root | StateBuilder.To<any>): StateBuilder
    update(params: ((old: StateTransformer.Params<T>) => StateTransformer.Params<T> | void) | StateTransformer.Params<T>, builder?: StateBuilder.Root | StateBuilder.To<any>): StateBuilder {
        if (!this.state) throw new Error(`To use update() from StateObjectSelector, 'state' must be defined.`);
        if (!builder) builder = this.state.build();
        (builder || this.state.build()).to(this).update(params);
        return builder;
    }

    /** Checks if the object exists. If not throw an error. */
    checkValid() {
        if (!this.state) {
            throw new Error('Unassigned State.');
        }
        const cell = this.cell;
        if (!cell) {
            throw new Error(`Not created at all. Did you await/then the corresponding state update?`);
        }
        if (cell.status === 'ok') return true;
        if (cell.status === 'error') throw new Error(cell.errorText);
        if (cell.obj === StateObject.Null) throw new Error('The object is Null.');
        throw new Error(`Unresolved. Did you await/then the corresponding state update?`);
    }

    get isOk() {
        const cell = this.cell;
        return cell && cell.status === 'ok' && cell.obj !== StateObject.Null;
    }

    constructor(public ref: StateTransform.Ref, public state?: State) {
    }
}

export namespace StateObjectSelector {
    export type Obj<S extends StateObjectSelector> = S extends StateObjectSelector<infer A> ? A : never
    export type Transformer<S extends StateObjectSelector> = S extends StateObjectSelector<any, infer T> ? T : never
}

export type StateObjectRef<S extends StateObject = StateObject> = StateObjectSelector<S> | StateObjectCell<S> | StateTransform.Ref

export namespace StateObjectRef {
    export function resolveRef<S extends StateObject>(ref?: StateObjectRef<S>): StateTransform.Ref | undefined {
        if (!ref) return;
        if (typeof ref === 'string') return ref;
        if (StateObjectCell.is(ref)) return ref.transform.ref;
        return ref.cell?.transform.ref;
    }

    export function resolve<S extends StateObject>(state: State, ref?: StateObjectRef<S>): StateObjectCell<S> | undefined {
        if (!ref) return;
        if (StateObjectCell.is(ref)) return ref;
        if (typeof ref === 'string') return state.cells.get(ref) as StateObjectCell<S> | undefined;
        return ref.cell;
    }

    export function resolveAndCheck<S extends StateObject>(state: State, ref?: StateObjectRef<S>): StateObjectCell<S> | undefined {
        const cell = resolve(state, ref);
        if (!cell || !cell.obj || cell.status !== 'ok') return;
        return cell;
    }
}