/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '../../../mol-model/structure';
import { PluginContext } from '../../../mol-plugin/context';
import { StateSelection, StateTransform, StateTransformer, StateObject, StateObjectCell } from '../../../mol-state';
import { StateTransforms } from '../../transforms';
import { PluginCommands } from '../../../mol-plugin/commands';
import { arraySetAdd } from '../../../mol-util/array';
import { PluginStateObject } from '../../objects';
import { StatefulPluginComponent } from '../../component';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { MeasurementRepresentationCommonTextParams } from '../../../mol-repr/shape/loci/common';

export { StructureMeasurementManager };

export const MeasurementGroupTag = 'measurement-group';

export type StructureMeasurementCell = StateObjectCell<PluginStateObject.Shape.Representation3D, StateTransform<StateTransformer<PluginStateObject.Molecule.Structure.Selections, PluginStateObject.Shape.Representation3D, any>>>

export const StructureMeasurementParams = {
    distanceUnitLabel: PD.Text('\u212B', { isEssential: true }),
    textColor: MeasurementRepresentationCommonTextParams.textColor
};
const DefaultStructureMeasurementOptions = PD.getDefaultValues(StructureMeasurementParams);
export type StructureMeasurementOptions = PD.ValuesFor<typeof StructureMeasurementParams>

export interface StructureMeasurementManagerState {
    labels: StructureMeasurementCell[],
    distances: StructureMeasurementCell[],
    angles: StructureMeasurementCell[],
    dihedrals: StructureMeasurementCell[],
    orientations: StructureMeasurementCell[],
    options: StructureMeasurementOptions
}

type StructureMeasurementManagerAddOptions = {
    customText?: string,
    selectionTags?: string | string[],
    reprTags?: string | string[]
}

class StructureMeasurementManager extends StatefulPluginComponent<StructureMeasurementManagerState>  {
    readonly behaviors = {
        state: this.ev.behavior(this.state)
    };

    private stateUpdated() {
        this.behaviors.state.next(this.state);
    }

    private getGroup() {
        const state = this.plugin.state.data;
        const groupRef = StateSelection.findTagInSubtree(state.tree, StateTransform.RootRef, MeasurementGroupTag);
        const builder = this.plugin.state.data.build();

        if (groupRef) return builder.to(groupRef);
        return builder.toRoot().group(StateTransforms.Misc.CreateGroup, { label: `Measurements` }, { tags: MeasurementGroupTag });
    }

    async setOptions(options: StructureMeasurementOptions) {
        if (this.updateState({ options })) this.stateUpdated();

        const update = this.plugin.state.data.build();
        for (const cell of this.state.distances) {
            update.to(cell).update((old: any) => {
                old.unitLabel = options.distanceUnitLabel;
                old.textColor = options.textColor;
            });
        }
        for (const cell of this.state.labels) {
            update.to(cell).update((old: any) => { old.textColor = options.textColor; });
        }
        for (const cell of this.state.angles) {
            update.to(cell).update((old: any) => { old.textColor = options.textColor; });
        }
        for (const cell of this.state.dihedrals) {
            update.to(cell).update((old: any) => { old.textColor = options.textColor; });
        }

        if (update.editInfo.count === 0) return;

        await PluginCommands.State.Update(this.plugin, { state: this.plugin.state.data, tree: update, options: { doNotLogTiming: true } });
    }

    async addDistance(a: StructureElement.Loci, b: StructureElement.Loci, options?: StructureMeasurementManagerAddOptions) {
        const cellA = this.plugin.helpers.substructureParent.get(a.structure);
        const cellB = this.plugin.helpers.substructureParent.get(b.structure);

        if (!cellA || !cellB) return;

        const dependsOn = [cellA.transform.ref];
        arraySetAdd(dependsOn, cellB.transform.ref);

        const update = this.getGroup();
        update
            .apply(StateTransforms.Model.MultiStructureSelectionFromExpression, {
                selections: [
                    { key: 'a', groupId: 'a', ref: cellA.transform.ref, expression: StructureElement.Loci.toExpression(a) },
                    { key: 'b', groupId: 'b', ref: cellB.transform.ref, expression: StructureElement.Loci.toExpression(b) }
                ],
                isTransitive: true,
                label: 'Distance'
            }, { dependsOn, tags: options?.selectionTags })
            .apply(StateTransforms.Representation.StructureSelectionsDistance3D, {
                customText: options?.customText || '',
                unitLabel: this.state.options.distanceUnitLabel,
                textColor: this.state.options.textColor
            }, { tags: options?.reprTags });

        const state = this.plugin.state.data;
        await PluginCommands.State.Update(this.plugin, { state, tree: update, options: { doNotLogTiming: true } });
    }

    async addAngle(a: StructureElement.Loci, b: StructureElement.Loci, c: StructureElement.Loci, options?: StructureMeasurementManagerAddOptions) {
        const cellA = this.plugin.helpers.substructureParent.get(a.structure);
        const cellB = this.plugin.helpers.substructureParent.get(b.structure);
        const cellC = this.plugin.helpers.substructureParent.get(c.structure);

        if (!cellA || !cellB || !cellC) return;

        const dependsOn = [cellA.transform.ref];
        arraySetAdd(dependsOn, cellB.transform.ref);
        arraySetAdd(dependsOn, cellC.transform.ref);

        const update = this.getGroup();
        update
            .apply(StateTransforms.Model.MultiStructureSelectionFromExpression, {
                selections: [
                    { key: 'a', ref: cellA.transform.ref, expression: StructureElement.Loci.toExpression(a) },
                    { key: 'b', ref: cellB.transform.ref, expression: StructureElement.Loci.toExpression(b) },
                    { key: 'c', ref: cellC.transform.ref, expression: StructureElement.Loci.toExpression(c) }
                ],
                isTransitive: true,
                label: 'Angle'
            }, { dependsOn, tags: options?.selectionTags })
            .apply(StateTransforms.Representation.StructureSelectionsAngle3D, {
                customText: options?.customText || '',
                textColor: this.state.options.textColor
            }, { tags: options?.reprTags });

        const state = this.plugin.state.data;
        await PluginCommands.State.Update(this.plugin, { state, tree: update, options: { doNotLogTiming: true } });
    }

    async addDihedral(a: StructureElement.Loci, b: StructureElement.Loci, c: StructureElement.Loci, d: StructureElement.Loci, options?: StructureMeasurementManagerAddOptions) {
        const cellA = this.plugin.helpers.substructureParent.get(a.structure);
        const cellB = this.plugin.helpers.substructureParent.get(b.structure);
        const cellC = this.plugin.helpers.substructureParent.get(c.structure);
        const cellD = this.plugin.helpers.substructureParent.get(d.structure);

        if (!cellA || !cellB || !cellC || !cellD) return;

        const dependsOn = [cellA.transform.ref];
        arraySetAdd(dependsOn, cellB.transform.ref);
        arraySetAdd(dependsOn, cellC.transform.ref);
        arraySetAdd(dependsOn, cellD.transform.ref);

        const update = this.getGroup();
        update
            .apply(StateTransforms.Model.MultiStructureSelectionFromExpression, {
                selections: [
                    { key: 'a', ref: cellA.transform.ref, expression: StructureElement.Loci.toExpression(a) },
                    { key: 'b', ref: cellB.transform.ref, expression: StructureElement.Loci.toExpression(b) },
                    { key: 'c', ref: cellC.transform.ref, expression: StructureElement.Loci.toExpression(c) },
                    { key: 'd', ref: cellD.transform.ref, expression: StructureElement.Loci.toExpression(d) }
                ],
                isTransitive: true,
                label: 'Dihedral'
            }, { dependsOn, tags: options?.selectionTags })
            .apply(StateTransforms.Representation.StructureSelectionsDihedral3D, {
                customText: options?.customText || '',
                textColor: this.state.options.textColor
            }, { tags: options?.reprTags });

        const state = this.plugin.state.data;
        await PluginCommands.State.Update(this.plugin, { state, tree: update, options: { doNotLogTiming: true } });
    }

    async addLabel(a: StructureElement.Loci, options?: Omit<StructureMeasurementManagerAddOptions, 'customText'>) {
        const cellA = this.plugin.helpers.substructureParent.get(a.structure);

        if (!cellA) return;

        const dependsOn = [cellA.transform.ref];

        const update = this.getGroup();
        update
            .apply(StateTransforms.Model.MultiStructureSelectionFromExpression, {
                selections: [
                    { key: 'a', ref: cellA.transform.ref, expression: StructureElement.Loci.toExpression(a) },
                ],
                isTransitive: true,
                label: 'Label'
            }, { dependsOn, tags: options?.selectionTags })
            .apply(StateTransforms.Representation.StructureSelectionsLabel3D, {
                textColor: this.state.options.textColor
            }, { tags: options?.reprTags });

        const state = this.plugin.state.data;
        await PluginCommands.State.Update(this.plugin, { state, tree: update, options: { doNotLogTiming: true } });
    }

    async addOrientation(a: StructureElement.Loci) {
        const cellA = this.plugin.helpers.substructureParent.get(a.structure);

        if (!cellA) return;

        const dependsOn = [cellA.transform.ref];

        const update = this.getGroup();
        update
            .apply(StateTransforms.Model.MultiStructureSelectionFromExpression, {
                selections: [
                    { key: 'a', ref: cellA.transform.ref, expression: StructureElement.Loci.toExpression(a) },
                ],
                isTransitive: true,
                label: 'Orientation'
            }, { dependsOn })
            .apply(StateTransforms.Representation.StructureSelectionsOrientation3D);

        const state = this.plugin.state.data;
        await PluginCommands.State.Update(this.plugin, { state, tree: update, options: { doNotLogTiming: true } });
    }

    private _empty: any[] = [];
    private getTransforms<T extends StateTransformer<A, B, any>, A extends PluginStateObject.Molecule.Structure.Selections, B extends StateObject>(transformer: T) {
        const state = this.plugin.state.data;
        const groupRef = StateSelection.findTagInSubtree(state.tree, StateTransform.RootRef, MeasurementGroupTag);
        const ret = groupRef ? state.select(StateSelection.Generators.ofTransformer(transformer, groupRef)) : this._empty;
        if (ret.length === 0) return this._empty;
        return ret;
    }

    private sync() {
        const updated = this.updateState({
            labels: this.getTransforms(StateTransforms.Representation.StructureSelectionsLabel3D),
            distances: this.getTransforms(StateTransforms.Representation.StructureSelectionsDistance3D),
            angles: this.getTransforms(StateTransforms.Representation.StructureSelectionsAngle3D),
            dihedrals: this.getTransforms(StateTransforms.Representation.StructureSelectionsDihedral3D),
            orientations: this.getTransforms(StateTransforms.Representation.StructureSelectionsOrientation3D)
        });
        if (updated) this.stateUpdated();
    }

    constructor(private plugin: PluginContext) {
        super({ labels: [], distances: [], angles: [], dihedrals: [], orientations: [], options: DefaultStructureMeasurementOptions });

        plugin.state.data.events.changed.subscribe(e => {
            if (e.inTransaction || plugin.behaviors.state.isAnimating.value) return;
            this.sync();
        });

        plugin.behaviors.state.isAnimating.subscribe(isAnimating => {
            if (!isAnimating && !plugin.behaviors.state.isUpdating.value) this.sync();
        });
    }
}