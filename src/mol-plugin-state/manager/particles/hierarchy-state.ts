/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject as SO } from '../../objects';
import type { StateTransforms } from '../../transforms';
import { StateObject, StateTransform, State, StateObjectCell, StateTree, StateTransformer } from '../../../mol-state';

export function buildParticleHierarchy(state: State, previous?: ParticleHierarchy) {
    const build = BuildState(state, previous || ParticleHierarchy());
    doPreOrder(state.tree, build);
    if (previous) previous.refs.forEach(isRemoved, build);
    return { hierarchy: build.hierarchy, added: build.added, changed: build.changed };
}

export interface ParticleHierarchy {
    lists: ParticleListRef[],
    refs: Map<StateTransform.Ref, ParticleHierarchyRef>
}

export function ParticleHierarchy(): ParticleHierarchy {
    return { lists: [], refs: new Map() };
}

interface RefBase<K extends string = string, O extends StateObject = StateObject, T extends StateTransformer = StateTransformer> {
    kind: K,
    cell: StateObjectCell<O, StateTransform<T>>,
    version: StateTransform['version']
}

export type ParticleHierarchyRef = ParticleListRef | ParticleRepresentationRef

export interface ParticleListRef extends RefBase<'particle-list', SO.Particle.List> {
    representations: ParticleRepresentationRef[]
}

function ParticleListRef(cell: StateObjectCell<SO.Particle.List>): ParticleListRef {
    return { kind: 'particle-list', cell, version: cell.transform.version, representations: [] };
}

export interface ParticleRepresentationRef extends RefBase<'particle-representation', SO.Particle.Representation3D, StateTransforms['Particles']['ParticlesRepresentation3D']> {
    list: ParticleListRef
}

function ParticleRepresentationRef(cell: StateObjectCell<SO.Particle.Representation3D>, list: ParticleListRef): ParticleRepresentationRef {
    return { kind: 'particle-representation', cell, version: cell.transform.version, list };
}

interface BuildState {
    state: State,
    oldHierarchy: ParticleHierarchy,

    hierarchy: ParticleHierarchy,

    currentList?: ParticleListRef,

    changed: boolean,
    added: Set<StateTransform.Ref>
}

function BuildState(state: State, oldHierarchy: ParticleHierarchy): BuildState {
    return { state, oldHierarchy, hierarchy: ParticleHierarchy(), changed: false, added: new Set() };
}

function createOrUpdateRefList<R extends ParticleHierarchyRef, C extends any[]>(state: BuildState, cell: StateObjectCell, list: R[], ctor: (...args: C) => R, ...args: C) {
    const ref: R = ctor(...args);
    list.push(ref);
    state.hierarchy.refs.set(cell.transform.ref, ref);
    const old = state.oldHierarchy.refs.get(cell.transform.ref);
    if (old) {
        if (old.version !== cell.transform.version) state.changed = true;
    } else {
        state.added.add(ref.cell.transform.ref);
        state.changed = true;
    }
    return ref;
}

type TestCell = (cell: StateObjectCell, state: BuildState) => boolean
type ApplyRef = (state: BuildState, cell: StateObjectCell) => boolean | void
type LeaveRef = (state: BuildState) => any

function isTypeRoot(t: StateObject.Ctor, target: (state: BuildState) => any): TestCell {
    return (cell, state) => !target(state) && t.is(cell.obj);
}

function noop() { }

const Mapping: [TestCell, ApplyRef, LeaveRef][] = [
    [isTypeRoot(SO.Particle.List, t => t.currentList), (state, cell) => {
        state.currentList = createOrUpdateRefList(state, cell, state.hierarchy.lists, ParticleListRef, cell);
    }, state => state.currentList = void 0],

    [(cell, state) => {
        return !cell.state.isGhost && !!state.currentList && SO.Particle.Representation3D.is(cell.obj);
    }, (state, cell) => {
        if (state.currentList) {
            createOrUpdateRefList(state, cell, state.currentList.representations, ParticleRepresentationRef, cell, state.currentList);
        }
        return false;
    }, noop]
];

function isValidCell(cell?: StateObjectCell): cell is StateObjectCell {
    if (!cell || !cell?.parent || !cell.parent.cells.has(cell.transform.ref)) return false;
    const { obj } = cell;
    if (!obj || obj === StateObject.Null || (cell.status !== 'ok' && cell.status !== 'error')) return false;
    return true;
}

function isRemoved(this: BuildState, ref: ParticleHierarchyRef) {
    const { cell } = ref;
    if (isValidCell(cell)) return;
    this.changed = true;
}

type VisitorCtx = { tree: StateTree, state: BuildState };

function _preOrderFunc(this: VisitorCtx, c: StateTransform.Ref | undefined) { _doPreOrder(this, this.tree.transforms.get(c!)!); }
function _doPreOrder(ctx: VisitorCtx, root: StateTransform) {
    const { state } = ctx;
    const cell = state.state.cells.get(root.ref);
    if (!isValidCell(cell)) return;

    let onLeave: undefined | ((state: BuildState) => any) = void 0;
    let end = false;
    for (const [test, f, l] of Mapping) {
        if (test(cell, state)) {
            const cont = f(state, cell);
            if (cont === false) {
                end = true;
                break;
            }
            onLeave = l;
            break;
        }
    }

    if (end) return;

    const children = ctx.tree.children.get(root.ref);
    if (children && children.size) {
        children.forEach(_preOrderFunc, ctx);
    }

    if (onLeave) onLeave(state);
}

function doPreOrder(tree: StateTree, state: BuildState): BuildState {
    const ctx: VisitorCtx = { tree, state };
    _doPreOrder(ctx, tree.root);
    return ctx.state;
}
