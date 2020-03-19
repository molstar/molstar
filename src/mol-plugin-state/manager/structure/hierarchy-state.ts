/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject as SO } from '../../objects';
import { StateObject, StateTransform, State, StateObjectCell, StateTree, StateTransformer } from '../../../mol-state';
import { StructureBuilderTags } from '../../builder/structure';
import { StructureRepresentationBuilderTags } from '../../builder/structure/representation';
import { StructureRepresentationInteractionTags } from '../../../mol-plugin/behavior/dynamic/selection/structure-representation-interaction';
import { StateTransforms } from '../../transforms';
import { VolumeStreaming } from '../../../mol-plugin/behavior/dynamic/volume-streaming/behavior';
import { CreateVolumeStreamingBehavior } from '../../../mol-plugin/behavior/dynamic/volume-streaming/transformers';

export function buildStructureHierarchy(state: State, previous?: StructureHierarchy) {
    const build = BuildState(state, previous || StructureHierarchy());
    doPreOrder(state.tree, build);
    if (previous) previous.refs.forEach(isRemoved, build);
    return { hierarchy: build.hierarchy, added: build.added, updated: build.updated, removed: build.removed };
}

export interface StructureHierarchy {
    trajectories: TrajectoryRef[],
    models: ModelRef[],
    structures: StructureRef[],
    refs: Map<StateTransform.Ref, HierarchyRef>
}

export function StructureHierarchy(): StructureHierarchy {
    return { trajectories: [], models: [], structures: [], refs: new Map() }
}

interface RefBase<K extends string = string, O extends StateObject = StateObject, T extends StateTransformer = StateTransformer> {
    kind: K,
    cell: StateObjectCell<O, StateTransform<T>>,
    version: StateTransform['version']
}

export type HierarchyRef =
    | TrajectoryRef
    | ModelRef | ModelPropertiesRef | ModelUnitcellRef
    | StructureRef | StructurePropertiesRef | StructureVolumeStreamingRef | StructureComponentRef | StructureRepresentationRef
    | GenericRepresentationRef

export interface TrajectoryRef extends RefBase<'trajectory', SO.Molecule.Trajectory> {
    models: ModelRef[]
}

function TrajectoryRef(cell: StateObjectCell<SO.Molecule.Trajectory>): TrajectoryRef {
    return { kind: 'trajectory', cell, version: cell.transform.version, models: [] };
}

export interface ModelRef extends RefBase<'model', SO.Molecule.Model> {
    trajectory?: TrajectoryRef,
    properties?: ModelPropertiesRef,
    structures: StructureRef[],
    genericRepresentations?: GenericRepresentationRef[],
    unitcell?: ModelUnitcellRef,
    /** to support decorators */
    childRoot: StateObjectCell<SO.Molecule.Model>
}

function ModelRef(cell: StateObjectCell<SO.Molecule.Model>, trajectory?: TrajectoryRef): ModelRef {
    return { kind: 'model', cell, version: cell.transform.version, trajectory, structures: [], childRoot: cell };
}

export interface ModelPropertiesRef extends RefBase<'model-properties', SO.Molecule.Model, StateTransforms['Model']['CustomModelProperties']> {
    model: ModelRef
}

function ModelPropertiesRef(cell: StateObjectCell<SO.Molecule.Model>, model: ModelRef): ModelPropertiesRef {
    return { kind: 'model-properties', cell, version: cell.transform.version, model };
}

export interface ModelUnitcellRef extends RefBase<'model-unitcell', SO.Shape.Representation3D, StateTransforms['Representation']['ModelUnitcell3D']> {
    model: ModelRef
}

function ModelUnitcellRef(cell: StateObjectCell<SO.Shape.Representation3D>, model: ModelRef): ModelUnitcellRef {
    return { kind: 'model-unitcell', cell, version: cell.transform.version, model };
}

export interface StructureRef extends RefBase<'structure', SO.Molecule.Structure> {
    model?: ModelRef,
    properties?: StructurePropertiesRef,
    components: StructureComponentRef[],
    currentFocus?: {
        focus?: StructureComponentRef,
        surroundings?: StructureComponentRef,
    },
    genericRepresentations?: GenericRepresentationRef[],
    volumeStreaming?: StructureVolumeStreamingRef,
    /** to support decorators */
    childRoot: StateObjectCell<SO.Molecule.Structure>
}

function StructureRef(cell: StateObjectCell<SO.Molecule.Structure>, model?: ModelRef): StructureRef {
    return { kind: 'structure', cell, version: cell.transform.version, model, components: [], childRoot: cell };
}

export interface StructurePropertiesRef extends RefBase<'structure-properties', SO.Molecule.Structure, StateTransforms['Model']['CustomStructureProperties']> {
    structure: StructureRef
}

function StructurePropertiesRef(cell: StateObjectCell<SO.Molecule.Structure>, structure: StructureRef): StructurePropertiesRef {
    return { kind: 'structure-properties', cell, version: cell.transform.version, structure };
}

export interface StructureVolumeStreamingRef extends RefBase<'structure-volume-streaming', VolumeStreaming, CreateVolumeStreamingBehavior> {
    structure: StructureRef
}

function StructureVolumeStreamingRef(cell: StateObjectCell<VolumeStreaming>, structure: StructureRef): StructureVolumeStreamingRef {
    return { kind: 'structure-volume-streaming', cell, version: cell.transform.version, structure };
}

export interface StructureComponentRef extends RefBase<'structure-component', SO.Molecule.Structure, StateTransforms['Model']['StructureComponent']> {
    structure: StructureRef,
    key?: string,
    representations: StructureRepresentationRef[],
    genericRepresentations?: GenericRepresentationRef[]
}

function componentKey(cell: StateObjectCell<SO.Molecule.Structure>) {
    if (!cell.transform.tags) return;
    return [...cell.transform.tags].sort().join();
}

function StructureComponentRef(cell: StateObjectCell<SO.Molecule.Structure>, structure: StructureRef): StructureComponentRef {
    return { kind: 'structure-component', cell, version: cell.transform.version, structure, key: componentKey(cell), representations: [] };
}

export interface StructureRepresentationRef extends RefBase<'structure-representation', SO.Molecule.Structure.Representation3D, StateTransforms['Representation']['StructureRepresentation3D']> {
    component: StructureComponentRef
}

function StructureRepresentationRef(cell: StateObjectCell<SO.Molecule.Structure.Representation3D>, component: StructureComponentRef): StructureRepresentationRef {
    return { kind: 'structure-representation', cell, version: cell.transform.version, component };
}

export interface GenericRepresentationRef extends RefBase<'generic-representation', SO.Any> {
    parent: HierarchyRef
}

function GenericRepresentationRef(cell: StateObjectCell<SO.Molecule.Structure.Representation3D>, parent: HierarchyRef): GenericRepresentationRef {
    return { kind: 'generic-representation', cell, version: cell.transform.version, parent };
}

interface BuildState {
    state: State,
    oldHierarchy: StructureHierarchy,

    hierarchy: StructureHierarchy,

    currentTrajectory?: TrajectoryRef,
    currentModel?: ModelRef,
    currentStructure?: StructureRef,
    currentComponent?: StructureComponentRef,

    updated: HierarchyRef[],
    added: HierarchyRef[],
    removed: HierarchyRef[]
}

function BuildState(state: State, oldHierarchy: StructureHierarchy): BuildState {
    return { state, oldHierarchy, hierarchy: StructureHierarchy(), updated: [], added: [], removed: [] };
}

function createOrUpdateRefList<R extends HierarchyRef, C extends any[]>(state: BuildState, cell: StateObjectCell, list: R[], ctor: (...args: C) => R, ...args: C) {
    const ref: R = ctor(...args);
    list.push(ref);
    state.hierarchy.refs.set(cell.transform.ref, ref);
    const old = state.oldHierarchy.refs.get(cell.transform.ref);
    if (old) {
        if (old.version !== cell.transform.version) state.updated.push(ref);
    } else {
        state.added.push(ref);
    }
    return ref;
}

function createOrUpdateRef<R extends HierarchyRef, C extends any[]>(state: BuildState, cell: StateObjectCell, ctor: (...args: C) => R, ...args: C) {
    const ref: R = ctor(...args);
    state.hierarchy.refs.set(cell.transform.ref, ref);
    const old = state.oldHierarchy.refs.get(cell.transform.ref);
    if (old) {
        if (old.version !== cell.transform.version) state.updated.push(ref);
    } else {
        state.added.push(ref);
    }
    return ref;
}

const tagMap: [string, (state: BuildState, cell: StateObjectCell) => boolean | void, (state: BuildState) => any][] = [
    [StructureBuilderTags.Trajectory, (state, cell) => {
        state.currentTrajectory = createOrUpdateRefList(state, cell, state.hierarchy.trajectories, TrajectoryRef, cell);
    }, state => state.currentTrajectory = void 0],
    [StructureBuilderTags.Model, (state, cell) => {
        if (state.currentTrajectory) {
            state.currentModel = createOrUpdateRefList(state, cell, state.currentTrajectory.models, ModelRef, cell, state.currentTrajectory);
        } else {
            state.currentModel = createOrUpdateRef(state, cell, ModelRef, cell);
        }
        state.hierarchy.models.push(state.currentModel);
    }, state => state.currentModel = void 0],
    [StructureBuilderTags.ModelProperties, (state, cell) => {
        if (!state.currentModel) return false;
        state.currentModel.properties = createOrUpdateRef(state, cell, ModelPropertiesRef, cell, state.currentModel);
    }, state => { }],
    [StructureBuilderTags.ModelUnitcell, (state, cell) => {
        if (!state.currentModel) return false;
        state.currentModel.unitcell = createOrUpdateRef(state, cell, ModelUnitcellRef, cell, state.currentModel);
    }, state => { }],
    [StructureBuilderTags.ModelGenericRepresentation, (state, cell) => {
        if (!state.currentModel) return false;
        if (!state.currentModel.genericRepresentations) state.currentModel.genericRepresentations = []
        createOrUpdateRefList(state, cell, state.currentModel.genericRepresentations, GenericRepresentationRef, cell, state.currentModel);
    }, state => { }],
    [StructureBuilderTags.Structure, (state, cell) => {
        if (state.currentModel) {
            state.currentStructure = createOrUpdateRefList(state, cell, state.currentModel.structures, StructureRef, cell, state.currentModel);
        } else {
            state.currentStructure = createOrUpdateRef(state, cell, StructureRef, cell);
        }
        state.hierarchy.structures.push(state.currentStructure);
    }, state => state.currentStructure = void 0],
    [StructureBuilderTags.StructureProperties, (state, cell) => {
        if (!state.currentStructure) return false;
        state.currentStructure.properties = createOrUpdateRef(state, cell, StructurePropertiesRef, cell, state.currentStructure);
    }, state => { }],
    [StructureBuilderTags.Component, (state, cell) => {
        if (!state.currentStructure) return false;
        state.currentComponent = createOrUpdateRefList(state, cell, state.currentStructure.components, StructureComponentRef, cell, state.currentStructure);
    }, state => state.currentComponent = void 0],
    [StructureRepresentationBuilderTags.Representation, (state, cell) => {
        if (!state.currentComponent) return false;
        createOrUpdateRefList(state, cell, state.currentComponent.representations, StructureRepresentationRef, cell, state.currentComponent);
    }, state => { }],
    [StructureRepresentationInteractionTags.ResidueSel, (state, cell) => {
        if (!state.currentStructure) return false;
        if (!state.currentStructure.currentFocus) state.currentStructure.currentFocus = { };
        state.currentStructure.currentFocus.focus = createOrUpdateRef(state, cell, StructureComponentRef, cell, state.currentStructure);
        state.currentComponent = state.currentStructure.currentFocus.focus;
    }, state => state.currentComponent = void 0],
    [StructureRepresentationInteractionTags.SurrSel, (state, cell) => {
        if (!state.currentStructure) return false;
        if (!state.currentStructure.currentFocus) state.currentStructure.currentFocus = { };
        state.currentStructure.currentFocus.surroundings = createOrUpdateRef(state, cell, StructureComponentRef, cell, state.currentStructure);
        state.currentComponent = state.currentStructure.currentFocus.surroundings;
    }, state => state.currentComponent = void 0]
]

function isValidCell(cell?: StateObjectCell): cell is StateObjectCell {
    if (!cell || !cell.parent.cells.has(cell.transform.ref)) return false;
    const { obj } = cell;
    if (!obj || obj === StateObject.Null || (cell.status !== 'ok' && cell.status !== 'error')) return false;
    return true;
}

function isRemoved(this: BuildState, ref: HierarchyRef) {
    const { cell } = ref;
    if (isValidCell(cell)) return;
    this.removed.push(ref);
}

type VisitorCtx = { tree: StateTree, state: BuildState };

function _preOrderFunc(this: VisitorCtx, c: StateTransform.Ref | undefined) { _doPreOrder(this, this.tree.transforms.get(c!)!); }
function _doPreOrder(ctx: VisitorCtx, root: StateTransform) {
    const { state } = ctx;
    const cell = state.state.cells.get(root.ref);
    if (!isValidCell(cell)) return;

    let onLeave: undefined | ((state: BuildState) => any) = void 0;
    for (const [t, f, l] of tagMap) {
        if (StateObject.hasTag(cell.obj!, t)) {
            const stop = f(state, cell);
            if (stop === false) {
                return;
            }
            onLeave = l;
            break;
        }
    }

    if (cell.transform.isDecorator) {
        if (state.currentModel && SO.Molecule.Model.is(cell.obj)) {
            state.currentModel.childRoot = cell;
        } else if (state.currentStructure && !state.currentComponent && SO.Molecule.Structure.is(cell.obj)) {
            state.currentStructure.childRoot = cell;
        }
    }

    if (!onLeave && !cell.state.isGhost && state.currentComponent && SO.Molecule.Structure.Representation3D.is(cell.obj)) {
        createOrUpdateRefList(state, cell, state.currentComponent.representations, StructureRepresentationRef, cell, state.currentComponent);
    } else if (!cell.state.isGhost && SO.isRepresentation3D(cell.obj)) {
        const genericTarget = state.currentComponent || state.currentModel || state.currentStructure;
        if (genericTarget) {
            if (!genericTarget.genericRepresentations) genericTarget.genericRepresentations = [];
            genericTarget.genericRepresentations.push(createOrUpdateRef(state, cell, GenericRepresentationRef, cell, genericTarget));
        }
    } else if (state.currentStructure && VolumeStreaming.is(cell.obj)) {
        state.currentStructure.volumeStreaming = createOrUpdateRef(state, cell, StructureVolumeStreamingRef, cell, state.currentStructure);
        return;
    }

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