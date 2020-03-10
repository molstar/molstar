/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject as SO } from '../../objects';
import { StateObject, StateTransform, State, StateObjectCell, StateTree } from '../../../mol-state';
import { StructureBuilderTags } from '../../builder/structure';
import { RepresentationProviderTags } from '../../builder/structure/provider';

export function buildStructureHierarchy(state: State, previous?: StructureHierarchy) {
    const build = BuildState(state, previous || StructureHierarchy());
    StateTree.doPreOrder(state.tree, state.tree.root, build, visitCell);
    if (previous) previous.refs.forEach(isRemoved, build);
    return { hierarchy: build.hierarchy, added: build.added, updated: build.updated, removed: build.removed };
}

export interface StructureHierarchy {
    trajectories: TrajectoryRef[],
    refs: Map<StateTransform.Ref, HierarchyRef>
}

export function StructureHierarchy(): StructureHierarchy {
    return { trajectories: [], refs: new Map() }
}

interface RefBase<K extends string = string, T extends StateObject = StateObject> {
    kind: K,
    cell: StateObjectCell<T>,
    version: StateTransform['version']
}

export type HierarchyRef = 
    | TrajectoryRef
    | ModelRef | ModelPropertiesRef
    | StructureRef | StructurePropertiesRef | StructureComponentRef | StructureRepresentationRef

export interface TrajectoryRef extends RefBase<'trajectory', SO.Molecule.Trajectory> {
    models: ModelRef[]
}

function TrajectoryRef(cell: StateObjectCell<SO.Molecule.Trajectory>): TrajectoryRef { 
    return { kind: 'trajectory', cell, version: cell.transform.version, models: [] };
}

export interface ModelRef extends RefBase<'model', SO.Molecule.Model> {
    trajectory: TrajectoryRef,
    properties?: ModelPropertiesRef,
    structures: StructureRef[]
}

function ModelRef(cell: StateObjectCell<SO.Molecule.Model>, trajectory: TrajectoryRef): ModelRef {
    return { kind: 'model', cell, version: cell.transform.version, trajectory, structures: [] };
}

export interface ModelPropertiesRef extends RefBase<'model-properties', SO.Molecule.Model> {
    model: ModelRef
}

function ModelPropertiesRef(cell: StateObjectCell<SO.Molecule.Model>, model: ModelRef): ModelPropertiesRef {
    return { kind: 'model-properties', cell, version: cell.transform.version, model };
}

export interface StructureRef extends RefBase<'structure', SO.Molecule.Structure> {
    model: ModelRef,
    properties?: StructurePropertiesRef,
    components: StructureComponentRef[],
    // representations: StructureRepresentationRef[],
    currentFocus?: {
        focus?: StructureComponentRef,
        surroundings?: StructureComponentRef,
    },
    // volumeStreaming?: ....
}

function StructureRef(cell: StateObjectCell<SO.Molecule.Structure>, model: ModelRef): StructureRef {
    return { kind: 'structure', cell, version: cell.transform.version, model, components: [] };
}

export interface StructurePropertiesRef extends RefBase<'structure-properties', SO.Molecule.Structure> {
    structure: StructureRef
}

function StructurePropertiesRef(cell: StateObjectCell<SO.Molecule.Structure>, structure: StructureRef): StructurePropertiesRef {
    return { kind: 'structure-properties', cell, version: cell.transform.version, structure };
}

export interface StructureComponentRef extends RefBase<'structure-component', SO.Molecule.Structure> {
    structure: StructureRef,
    representations: StructureRepresentationRef[],
}

function StructureComponentRef(cell: StateObjectCell<SO.Molecule.Structure>, structure: StructureRef): StructureComponentRef {
    return { kind: 'structure-component', cell, version: cell.transform.version, structure, representations: [] };
}

export interface StructureRepresentationRef extends RefBase<'structure-representation', SO.Molecule.Structure.Representation3D> {
    component: StructureComponentRef
}

function StructureRepresentationRef(cell: StateObjectCell<SO.Molecule.Structure.Representation3D>, component: StructureComponentRef): StructureRepresentationRef {
    return { kind: 'structure-representation', cell, version: cell.transform.version, component };
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

function createOrUpdateRef<R extends HierarchyRef, C extends any[]>(state: BuildState, cell: StateObjectCell, old: R | undefined, ctor: (...args: C) => R, ...args: C) {
    const ref: R = ctor(...args);
    state.hierarchy.refs.set(cell.transform.ref, ref);
    if (old) {
        if (old.version !== cell.transform.version) state.updated.push(ref);
    } else {
        state.added.push(ref);
    }
    return ref;
}

const tagMap: [string, (state: BuildState, cell: StateObjectCell) => boolean | void][] = [
    [StructureBuilderTags.Trajectory, (state, cell) => {
        state.currentTrajectory = createOrUpdateRefList(state, cell, state.hierarchy.trajectories, TrajectoryRef, cell);
    }],
    [StructureBuilderTags.Model, (state, cell) => {
        if (!state.currentTrajectory) return false;
        state.currentModel = createOrUpdateRefList(state, cell, state.currentTrajectory.models, ModelRef, cell, state.currentTrajectory);
    }],
    [StructureBuilderTags.ModelProperties, (state, cell) => {
        if (!state.currentModel) return false;
        state.currentModel.properties = createOrUpdateRef(state, cell, state.currentModel.properties, ModelPropertiesRef, cell, state.currentModel);
    }],
    [StructureBuilderTags.Structure, (state, cell) => {
        if (!state.currentModel) return false;
        state.currentStructure = createOrUpdateRefList(state, cell, state.currentModel.structures, StructureRef, cell, state.currentModel);
    }],
    [StructureBuilderTags.StructureProperties, (state, cell) => {
        if (!state.currentStructure) return false;
        state.currentStructure.properties = createOrUpdateRef(state, cell, state.currentStructure.properties, StructurePropertiesRef, cell, state.currentStructure);
    }],
    [StructureBuilderTags.Component, (state, cell) => {
        if (!state.currentStructure) return false;
        state.currentComponent = createOrUpdateRefList(state, cell, state.currentStructure.components, StructureComponentRef, cell, state.currentStructure);
    }],
    [RepresentationProviderTags.Representation, (state, cell) => {
        if (!state.currentComponent) return false;
        createOrUpdateRefList(state, cell, state.currentComponent.representations, StructureRepresentationRef, cell, state.currentComponent);
    }]
]

function isValidCell(cell?: StateObjectCell): cell is StateObjectCell {
    if (!cell || !cell.parent.cells.has(cell.transform.ref)) return false;
    const { obj } = cell;
    if (!obj || obj === StateObject.Null || (cell.status !== 'ok' && cell.status !== 'error')) return false;
    return true;
}

function visitCell(t: StateTransform, tree: StateTree, state: BuildState): boolean {
    const cell = state.state.cells.get(t.ref);
    if (!isValidCell(cell)) return false;

    for (const [t, f] of tagMap) {
        if (StateObject.hasTag(cell.obj!, t)) {
            const stop = f(state, cell);
            if (stop === false) return false;
            return true;
        }
    }

    if (state.currentComponent && SO.Molecule.Structure.Representation3D.is(cell.obj)) {
        createOrUpdateRefList(state, cell, state.currentComponent.representations, StructureRepresentationRef, cell, state.currentComponent);
    }

    return true;
}

function isRemoved(this: BuildState, ref: HierarchyRef) {
    const { cell } = ref;
    if (isValidCell(cell)) return;
    this.removed.push(ref);
}
