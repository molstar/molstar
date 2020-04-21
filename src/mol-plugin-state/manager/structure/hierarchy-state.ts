/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject as SO } from '../../objects';
import { StateObject, StateTransform, State, StateObjectCell, StateTree, StateTransformer } from '../../../mol-state';
import { StateTransforms } from '../../transforms';
import { VolumeStreaming } from '../../../mol-plugin/behavior/dynamic/volume-streaming/behavior';
import { CreateVolumeStreamingBehavior } from '../../../mol-plugin/behavior/dynamic/volume-streaming/transformers';

export function buildStructureHierarchy(state: State, previous?: StructureHierarchy) {
    const build = BuildState(state, previous || StructureHierarchy());
    doPreOrder(state.tree, build);
    if (previous) previous.refs.forEach(isRemoved, build);
    return { hierarchy: build.hierarchy, added: build.added, changed: build.changed };
}

export interface StructureHierarchy {
    trajectories: TrajectoryRef[],
    models: ModelRef[],
    structures: StructureRef[],
    refs: Map<StateTransform.Ref, StructureHierarchyRef>
    // TODO: might be needed in the future
    // decorators: Map<StateTransform.Ref, StateTransform>,
}

export function StructureHierarchy(): StructureHierarchy {
    return { trajectories: [], models: [], structures: [], refs: new Map() };
}

interface RefBase<K extends string = string, O extends StateObject = StateObject, T extends StateTransformer = StateTransformer> {
    kind: K,
    cell: StateObjectCell<O, StateTransform<T>>,
    version: StateTransform['version']
}

export type StructureHierarchyRef =
    | TrajectoryRef
    | ModelRef | ModelPropertiesRef | ModelUnitcellRef
    | StructureRef | StructurePropertiesRef | StructureTransformRef | StructureVolumeStreamingRef | StructureComponentRef | StructureRepresentationRef
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
    unitcell?: ModelUnitcellRef
}

function ModelRef(cell: StateObjectCell<SO.Molecule.Model>, trajectory?: TrajectoryRef): ModelRef {
    return { kind: 'model', cell, version: cell.transform.version, trajectory, structures: [] };
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
    transform?: StructureTransformRef,
    components: StructureComponentRef[],
    genericRepresentations?: GenericRepresentationRef[],
    volumeStreaming?: StructureVolumeStreamingRef
}

function StructureRef(cell: StateObjectCell<SO.Molecule.Structure>, model?: ModelRef): StructureRef {
    return { kind: 'structure', cell, version: cell.transform.version, model, components: [] };
}

export interface StructurePropertiesRef extends RefBase<'structure-properties', SO.Molecule.Structure, StateTransforms['Model']['CustomStructureProperties']> {
    structure: StructureRef
}

function StructurePropertiesRef(cell: StateObjectCell<SO.Molecule.Structure>, structure: StructureRef): StructurePropertiesRef {
    return { kind: 'structure-properties', cell, version: cell.transform.version, structure };
}

export interface StructureTransformRef extends RefBase<'structure-transform', SO.Molecule.Structure, StateTransforms['Model']['TransformStructureConformation']> {
    structure: StructureRef
}

function StructureTransformRef(cell: StateObjectCell<SO.Molecule.Structure>, structure: StructureRef): StructureTransformRef {
    return { kind: 'structure-transform', cell, version: cell.transform.version, structure };
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
    if (!cell.transform.tags) return cell.transform.ref;
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
    parent: StructureHierarchyRef
}

function GenericRepresentationRef(cell: StateObjectCell<SO.Molecule.Structure.Representation3D>, parent: StructureHierarchyRef): GenericRepresentationRef {
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

    changed: boolean,
    added: Set<StateTransform.Ref>
}

function BuildState(state: State, oldHierarchy: StructureHierarchy): BuildState {
    return { state, oldHierarchy, hierarchy: StructureHierarchy(), changed: false, added: new Set() };
}

function createOrUpdateRefList<R extends StructureHierarchyRef, C extends any[]>(state: BuildState, cell: StateObjectCell, list: R[], ctor: (...args: C) => R, ...args: C) {
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

function createOrUpdateRef<R extends StructureHierarchyRef, C extends any[]>(state: BuildState, cell: StateObjectCell, ctor: (...args: C) => R, ...args: C) {
    const ref: R = ctor(...args);
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

function isType(t: StateObject.Ctor): TestCell {
    return (cell) => t.is(cell.obj);
}

function isTypeRoot(t: StateObject.Ctor, target: (state: BuildState) => any): TestCell {
    return (cell, state) => !target(state) && t.is(cell.obj);
}

function isTransformer(t: StateTransformer): TestCell {
    return cell => cell.transform.transformer === t;
}

function noop() { }

const Mapping: [TestCell, ApplyRef, LeaveRef][] = [
    // Trajectory
    [isType(SO.Molecule.Trajectory), (state, cell) => {
        state.currentTrajectory = createOrUpdateRefList(state, cell, state.hierarchy.trajectories, TrajectoryRef, cell);
    }, state => state.currentTrajectory = void 0],

    // Model
    [isTypeRoot(SO.Molecule.Model, s => s.currentModel), (state, cell) => {
        if (state.currentTrajectory) {
            state.currentModel = createOrUpdateRefList(state, cell, state.currentTrajectory.models, ModelRef, cell, state.currentTrajectory);
        } else {
            state.currentModel = createOrUpdateRef(state, cell, ModelRef, cell);
        }
        state.hierarchy.models.push(state.currentModel);
    }, state => state.currentModel = void 0],
    [isTransformer(StateTransforms.Model.CustomModelProperties), (state, cell) => {
        if (!state.currentModel) return false;
        state.currentModel.properties = createOrUpdateRef(state, cell, ModelPropertiesRef, cell, state.currentModel);
    }, noop],
    [isTransformer(StateTransforms.Representation.ModelUnitcell3D), (state, cell) => {
        if (!state.currentModel) return false;
        state.currentModel.unitcell = createOrUpdateRef(state, cell, ModelUnitcellRef, cell, state.currentModel);
    }, noop],

    // Structure
    [isTypeRoot(SO.Molecule.Structure, s => s.currentStructure), (state, cell) => {
        if (state.currentModel) {
            state.currentStructure = createOrUpdateRefList(state, cell, state.currentModel.structures, StructureRef, cell, state.currentModel);
        } else {
            state.currentStructure = createOrUpdateRef(state, cell, StructureRef, cell);
        }
        state.hierarchy.structures.push(state.currentStructure);
    }, state => state.currentStructure = void 0],
    [isTransformer(StateTransforms.Model.CustomStructureProperties), (state, cell) => {
        if (!state.currentStructure) return false;
        state.currentStructure.properties = createOrUpdateRef(state, cell, StructurePropertiesRef, cell, state.currentStructure);
    }, noop],
    [isTransformer(StateTransforms.Model.TransformStructureConformation), (state, cell) => {
        if (!state.currentStructure) return false;
        state.currentStructure.transform = createOrUpdateRef(state, cell, StructureTransformRef, cell, state.currentStructure);
    }, noop],

    // Volume Streaming
    [isType(VolumeStreaming), (state, cell) => {
        if (!state.currentStructure) return false;
        state.currentStructure.volumeStreaming = createOrUpdateRef(state, cell, StructureVolumeStreamingRef, cell, state.currentStructure);
        // Do not continue into VolumeStreaming subtree.
        return false;
    }, noop],

    // Component
    [(cell, state) => {
        if (state.currentComponent || !state.currentStructure || cell.transform.transformer.definition.isDecorator) return false;
        return SO.Molecule.Structure.is(cell.obj);
    }, (state, cell) => {
        if (state.currentStructure) {
            state.currentComponent = createOrUpdateRefList(state, cell, state.currentStructure.components, StructureComponentRef, cell, state.currentStructure);
        }
    }, state => state.currentComponent = void 0],

    // Component Representation
    [(cell, state) => {
        return !cell.state.isGhost && !!state.currentComponent && SO.Molecule.Structure.Representation3D.is(cell.obj);
    }, (state, cell) => {
        if (state.currentComponent) {
            createOrUpdateRefList(state, cell, state.currentComponent.representations, StructureRepresentationRef, cell, state.currentComponent);
        }

        // Nothing useful down the line
        return false;
    }, noop],

    // Generic Representation
    [cell => !cell.state.isGhost && SO.isRepresentation3D(cell.obj), (state, cell) => {
        const genericTarget = state.currentComponent || state.currentStructure || state.currentModel;
        if (genericTarget) {
            if (!genericTarget.genericRepresentations) genericTarget.genericRepresentations = [];
            createOrUpdateRefList(state, cell, genericTarget.genericRepresentations, GenericRepresentationRef, cell, genericTarget);
        }
    }, noop],
];

function isValidCell(cell?: StateObjectCell): cell is StateObjectCell {
    if (!cell || !cell?.parent || !cell.parent.cells.has(cell.transform.ref)) return false;
    const { obj } = cell;
    if (!obj || obj === StateObject.Null || (cell.status !== 'ok' && cell.status !== 'error')) return false;
    return true;
}

function isRemoved(this: BuildState, ref: StructureHierarchyRef) {
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