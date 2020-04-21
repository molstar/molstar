/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject as SO } from '../../objects';
import { StateObject, StateTransform, State, StateObjectCell, StateTree, StateTransformer } from '../../../mol-state';
import { StateTransforms } from '../../transforms';

export function buildVolumeHierarchy(state: State, previous?: VolumeHierarchy) {
    const build = BuildState(state, previous || VolumeHierarchy());
    doPreOrder(state.tree, build);
    if (previous) previous.refs.forEach(isRemoved, build);
    return { hierarchy: build.hierarchy, added: build.added, changed: build.changed };
}

export interface VolumeHierarchy {
    volumes: VolumeRef[],
    refs: Map<StateTransform.Ref, VolumeHierarchyRef>
    // TODO: might be needed in the future
    // decorators: Map<StateTransform.Ref, StateTransform>,
}

export function VolumeHierarchy(): VolumeHierarchy {
    return { volumes: [], refs: new Map() };
}

interface RefBase<K extends string = string, O extends StateObject = StateObject, T extends StateTransformer = StateTransformer> {
    kind: K,
    cell: StateObjectCell<O, StateTransform<T>>,
    version: StateTransform['version']
}

export type VolumeHierarchyRef = VolumeRef | VolumeRepresentationRef

export interface VolumeRef extends RefBase<'volume', SO.Volume.Data> {
    representations: VolumeRepresentationRef[]
}

function VolumeRef(cell: StateObjectCell<SO.Volume.Data>): VolumeRef {
    return { kind: 'volume', cell, version: cell.transform.version, representations: [] };
}

export interface VolumeRepresentationRef extends RefBase<'volume-representation', SO.Volume.Representation3D, StateTransforms['Representation']['VolumeRepresentation3D']> {
    volume: VolumeRef
}

function VolumeRepresentationRef(cell: StateObjectCell<SO.Volume.Representation3D>, volume: VolumeRef): VolumeRepresentationRef {
    return { kind: 'volume-representation', cell, version: cell.transform.version, volume };
}

interface BuildState {
    state: State,
    oldHierarchy: VolumeHierarchy,

    hierarchy: VolumeHierarchy,

    currentVolume?: VolumeRef,

    changed: boolean,
    added: Set<StateTransform.Ref>
}

function BuildState(state: State, oldHierarchy: VolumeHierarchy): BuildState {
    return { state, oldHierarchy, hierarchy: VolumeHierarchy(), changed: false, added: new Set() };
}

function createOrUpdateRefList<R extends VolumeHierarchyRef, C extends any[]>(state: BuildState, cell: StateObjectCell, list: R[], ctor: (...args: C) => R, ...args: C) {
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
    [isTypeRoot(SO.Volume.Data, t => t.currentVolume), (state, cell) => {
        state.currentVolume = createOrUpdateRefList(state, cell, state.hierarchy.volumes, VolumeRef, cell);
    }, state => state.currentVolume = void 0],

    [(cell, state) => {
        return !cell.state.isGhost && !!state.currentVolume && SO.Volume.Representation3D.is(cell.obj);
    }, (state, cell) => {
        if (state.currentVolume) {
            createOrUpdateRefList(state, cell, state.currentVolume.representations, VolumeRepresentationRef, cell, state.currentVolume);
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

function isRemoved(this: BuildState, ref: VolumeHierarchyRef) {
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

    // TODO: might be needed in the future
    // const { currentComponent, currentModel, currentStructure, currentTrajectory } = ctx.state;
    // const inTrackedSubtree = currentComponent || currentModel || currentStructure || currentTrajectory;

    // if (inTrackedSubtree && cell.transform.transformer.definition.isDecorator) {
    //     const ref = cell.transform.ref;
    //     const old = ctx.state.oldHierarchy.decorators.get(ref);
    //     if (old && old.version !== cell.transform.version) {
    //         ctx.state.changed = true;
    //     }
    //     ctx.state.hierarchy.decorators.set(cell.transform.ref, cell.transform);
    // }

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