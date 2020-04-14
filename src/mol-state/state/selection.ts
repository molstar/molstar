/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject, StateObjectCell } from '../object';
import { State } from '../state';
import { StateTree } from '../tree';
import { StateTransform } from '../transform';
import { StateTransformer } from '../transformer';

namespace StateSelection {
    export type Selector<C extends StateObjectCell = StateObjectCell> = Query<C> | Builder<C> | string | C;
    export type CellSeq<C extends StateObjectCell = StateObjectCell> = C[]
    export type Query<C extends StateObjectCell = StateObjectCell> = (state: State) => CellSeq<C>;

    export function select<C extends StateObjectCell>(s: Selector<C>, state: State) {
        return compile(s)(state);
    }

    export function compile<C extends StateObjectCell>(s: Selector<C>): Query<C> {
        const selector = s ? s : Generators.root;
        let query: Query;
        if (isBuilder(selector)) query = (selector as any).compile();
        else if (isObj(selector)) query = (Generators.byValue(selector) as any).compile();
        else if (isQuery(selector)) query = selector;
        else query = (Generators.byRef(selector as string) as any).compile();
        return query as Query<C>;
    }

    function isObj(arg: any): arg is StateObjectCell {
        return (arg as StateObjectCell).transform !== void 0 && (arg as StateObjectCell).status !== void 0;
    }

    function isBuilder(arg: any): arg is Builder {
        return arg.compile !== void 0;
    }

    function isQuery(arg: any): arg is Query {
        return typeof arg === 'function';
    }

    export interface Builder<C extends StateObjectCell = StateObjectCell> {
        flatMap<D extends StateObjectCell>(f: (n: C) => D[]): Builder<D>;
        mapObject<D extends StateObjectCell>(f: (n: C) => D): Builder<D>;
        unique(): Builder<C>;

        parent(): Builder<C>;
        first(): Builder<C>;
        filter(p: (n: C) => boolean): Builder<C>;
        withTag<D extends StateObjectCell = C>(tag: string): Builder<D>;
        withTransformer<T extends StateTransformer<any, StateObjectCell.Obj<C>, any>>(t: T): Builder<StateObjectCell<StateObjectCell.Obj<C>, StateTransform<T>>>;
        withStatus(s: StateObjectCell.Status): Builder<C>;
        subtree(): Builder;
        children(): Builder;
        ofType<T extends StateObject.Ctor>(t: T): Builder<StateObjectCell<StateObject.From<T>>>;
        ancestorOfType<T extends StateObject.Ctor>(t: T[]): Builder<StateObjectCell<StateObject.From<T>>>;
        rootOfType(t: StateObject.Ctor[]): Builder;

        select(state: State): CellSeq<C>
    }

    const BuilderPrototype: any = {
        select(state?: State) {
            return select(this, state || this.state);
        }
    };

    function registerModifier(name: string, f: Function) {
        BuilderPrototype[name] = function (this: any, ...args: any[]) { return f.call(void 0, this, ...args); };
    }

    function build<C extends StateObjectCell>(compile: () => Query<C>): Builder<C> {
        return Object.create(BuilderPrototype, { compile: { writable: false, configurable: false, value: compile } });
    }

    export namespace Generators {
        export const root = build(() => (state: State) => [state.cells.get(state.tree.root.ref)!]);

        export function byRef<T extends StateObject.Ctor>(...refs: StateTransform.Ref[]) {
            return build(() => (state: State) => {
                const ret: StateObjectCell<StateObject.From<T>>[] = [];
                for (const ref of refs) {
                    const n = state.cells.get(ref);
                    if (!n) continue;
                    ret.push(n as any);
                }
                return ret;
            });
        }

        export function byValue<T extends StateObjectCell>(...objects: T[]) { return build(() => (state: State) => objects); }

        export function rootsOfType<T extends StateObject.Ctor>(type: T, root: StateTransform.Ref = StateTransform.RootRef) {
            return build(() => state => {
                const ctx = { roots: [] as StateObjectCell<StateObject.From<T>>[], cells: state.cells, type: type.type };
                StateTree.doPreOrder(state.tree, state.tree.transforms.get(root), ctx, _findRootsOfType);
                return ctx.roots;
            });
        }

        export function ofType<T extends StateObject.Ctor>(type: T, root: StateTransform.Ref = StateTransform.RootRef) {
            return build(() => state => {
                const ctx = { ret: [] as StateObjectCell<StateObject.From<T>>[], cells: state.cells, type: type.type };
                StateTree.doPreOrder(state.tree, state.tree.transforms.get(root), ctx, _findOfType);
                return ctx.ret;
            });
        }

        export function ofTransformer<T extends StateTransformer<any, A, any>, A extends StateObject>(t: T, root: StateTransform.Ref = StateTransform.RootRef) {
            return build(() => state => {
                const ctx = { ret: [] as StateObjectCell<A, StateTransform<T>>[], cells: state.cells, t };
                StateTree.doPreOrder(state.tree, state.tree.transforms.get(root), ctx, _findOfTransformer);
                return ctx.ret;
            });
        }

        export function ofTransformerWithError<T extends StateTransformer<any, A, any>, A extends StateObject>(t: T, root: StateTransform.Ref = StateTransform.RootRef) {
            return build(() => state => {
                const ctx = { ret: [] as StateObjectCell<A, StateTransform<T>>[], cells: state.cells, t };
                StateTree.doPreOrder(state.tree, state.tree.transforms.get(root), ctx, _findOfTransformerWithError);
                return ctx.ret;
            });
        }

        function _findRootsOfType(n: StateTransform, _: any, s: { type: StateObject.Type, roots: StateObjectCell[], cells: State.Cells }) {
            const cell = s.cells.get(n.ref);
            if (cell && cell.obj && cell.obj.type === s.type) {
                s.roots.push(cell);
                return false;
            }
            return true;
        }

        function _findOfType(n: StateTransform, _: any, s: { type: StateObject.Type, ret: StateObjectCell[], cells: State.Cells }) {
            const cell = s.cells.get(n.ref);
            if (cell && cell.obj && cell.obj.type === s.type) {
                s.ret.push(cell);
            }
            return true;
        }

        function _findOfTransformer(n: StateTransform, _: any, s: { t: StateTransformer, ret: StateObjectCell[], cells: State.Cells }) {
            const cell = s.cells.get(n.ref);
            if (cell && cell.obj && cell.transform.transformer === s.t) {
                s.ret.push(cell);
            }
            return true;
        }

        function _findOfTransformerWithError(n: StateTransform, _: any, s: { t: StateTransformer, ret: StateObjectCell[], cells: State.Cells }) {
            const cell = s.cells.get(n.ref);
            if (cell && cell.status === 'error' && cell.transform.transformer === s.t) {
                s.ret.push(cell);
            }
            return true;
        }
    }

    registerModifier('flatMap', flatMap);
    export function flatMap(b: Selector, f: (obj: StateObjectCell, state: State) => CellSeq) {
        const q = compile(b);
        return build(() => (state: State) => {
            const ret: StateObjectCell[] = [];
            for (const n of q(state)) {
                for (const m of f(n, state)) {
                    ret.push(m);
                }
            }
            return ret;
        });
    }

    registerModifier('mapObject', mapObject);
    export function mapObject(b: Selector, f: (n: StateObjectCell, state: State) => StateObjectCell | undefined) {
        const q = compile(b);
        return build(() => (state: State) => {
            const ret: StateObjectCell[] = [];
            for (const n of q(state)) {
                const x = f(n, state);
                if (x) ret.push(x);
            }
            return ret;
        });
    }

    registerModifier('unique', unique);
    export function unique(b: Selector) {
        const q = compile(b);
        return build(() => (state: State) => {
            const set = new Set<string>();
            const ret: StateObjectCell[] = [];
            for (const n of q(state)) {
                if (!n) continue;
                if (!set.has(n.transform.ref)) {
                    set.add(n.transform.ref);
                    ret.push(n);
                }
            }
            return ret;
        });
    }

    registerModifier('first', first);
    export function first(b: Selector) {
        const q = compile(b);
        return build(() => (state: State) => {
            const r = q(state);
            return r.length ? [r[0]] : [];
        });
    }

    registerModifier('filter', filter);
    export function filter(b: Selector, p: (n: StateObjectCell) => boolean) { return flatMap(b, n => p(n) ? [n] : []); }

    registerModifier('withStatus', withStatus);
    export function withStatus(b: Selector, s: StateObjectCell.Status) { return filter(b, n => n.status === s); }

    registerModifier('withTag', withTag);
    export function withTag(b: Selector, tag: string) { return filter(b, n => !!n.transform.tags && n.transform.tags.indexOf(tag) >= 0); }

    registerModifier('subtree', subtree);
    export function subtree(b: Selector) {
        return flatMap(b, (n, s) => {
            const nodes = [] as string[];
            StateTree.doPreOrder(s.tree, s.tree.transforms.get(n.transform.ref), nodes, (x, _, ctx) => { ctx.push(x.ref); });
            return nodes.map(x => s.cells.get(x)!);
        });
    }

    registerModifier('children', children);
    export function children(b: Selector) {
        return flatMap(b, (n, s) => {
            const nodes: StateObjectCell[] = [];
            s.tree.children.get(n.transform.ref).forEach(c => nodes.push(s.cells.get(c!)!));
            return nodes;
        });
    }

    registerModifier('ofType', ofType);
    export function ofType(b: Selector, t: StateObject.Ctor) { return filter(b, n => n.obj ? n.obj.type === t.type : false); }

    registerModifier('ancestorOfType', ancestorOfType);
    export function ancestorOfType(b: Selector, types: StateObject.Ctor[]) { return unique(mapObject(b, (n, s) => findAncestorOfType(s.tree, s.cells, n.transform.ref, types))); }

    registerModifier('withTransformer', withTransformer);
    export function withTransformer(b: Selector, t: StateTransformer) { return filter(b, o => o.transform.transformer === t); }

    registerModifier('rootOfType', rootOfType);
    export function rootOfType(b: Selector, types: StateObject.Ctor[]) { return unique(mapObject(b, (n, s) => findRootOfType(s.tree, s.cells, n.transform.ref, types))); }

    registerModifier('parent', parent);
    export function parent(b: Selector) { return unique(mapObject(b, (n, s) => s.cells.get(s.tree.transforms.get(n.transform.ref)!.parent))); }

    export function findAncestorOfType<T extends StateObject.Ctor>(tree: StateTree, cells: State.Cells, root: StateTransform.Ref, types: T[]): StateObjectCell<StateObject.From<T>> | undefined {
        let current = tree.transforms.get(root)!, len = types.length;
        while (true) {
            current = tree.transforms.get(current.parent)!;
            const cell = cells.get(current.ref)!;
            if (!cell.obj) return void 0;
            const obj = cell.obj;
            for (let i = 0; i < len; i++) {
                if (obj.type === types[i].type) return cell as StateObjectCell<StateObject.From<T>>;
            }
            if (current.ref === StateTransform.RootRef) {
                return void 0;
            }
        }
    }

    export function findRootOfType(tree: StateTree, cells: State.Cells, root: StateTransform.Ref, types: StateObject.Ctor[]): StateObjectCell | undefined {
        let parent: StateObjectCell | undefined, _root = root;
        while (true) {
            const _parent = StateSelection.findAncestorOfType(tree, cells, _root, types);
            if (_parent) {
                parent = _parent;
                _root = _parent.transform.ref;
            } else {
                break;
            }
        }
        return parent;
    }

    export function findUniqueTagsInSubtree<K extends string = string>(tree: StateTree, root: StateTransform.Ref, tags: Set<K>): { [name in K]?: StateTransform.Ref } {
        return StateTree.doPreOrder(tree, tree.transforms.get(root), { refs: { }, tags }, _findUniqueTagsInSubtree).refs;
    }

    function _findUniqueTagsInSubtree(n: StateTransform, _: any, s: { refs: { [name: string]: StateTransform.Ref }, tags: Set<string> }) {
        if (n.tags) {
            for (const t of n.tags) {
                if (!s.tags.has(t)) continue;
                s.refs[t] = n.ref;
                break;
            }
        }
        return true;
    }

    export function findTagInSubtree(tree: StateTree, root: StateTransform.Ref, tag: string): StateTransform.Ref | undefined {
        return StateTree.doPreOrder(tree, tree.transforms.get(root), { ref: void 0, tag }, _findTagInSubtree).ref;
    }

    function _findTagInSubtree(n: StateTransform, _: any, s: { ref: string | undefined, tag: string }) {
        if (n.tags && n.tags.indexOf(s.tag) >= 0) {
            s.ref = n.ref;
            return false;
        }
        return true;
    }

    export function findWithAllTags<K extends string = string>(tree: StateTree, root: StateTransform.Ref, tags: Set<K>): StateTransform[] {
        return StateTree.doPreOrder(tree, tree.transforms.get(root), { refs: [], tags }, _findWithAllTags).refs;
    }

    function _findWithAllTags(n: StateTransform, _: any, s: { refs: StateTransform[], tags: Set<string> }) {
        if (n.tags) {
            const len = s.tags.size;
            let found = 0;
            for (const t of n.tags) {
                if (!s.tags.has(t)) continue;
                found++;

                if (found === len) {
                    s.refs.push(n);
                    break;
                }
            }
        } else if (s.tags.size === 0) {
            s.refs.push(n);
        }
    }

    export function tryFindDecorator<T extends StateTransformer>(state: State, root: StateTransform.Ref, transformer: T): StateObjectCell<StateTransformer.To<T>, StateTransform<T>> | undefined {
        const t = state.transforms.get(root);
        if (t.transformer === transformer) return state.cells.get(root)! as any;

        const children = state.tree.children.get(root);
        if (children.size !== 1) return;
        const first = children.first();
        if (state.transforms.get(first).transformer.definition.isDecorator) return tryFindDecorator(state, first, transformer);
    }
}

export { StateSelection };