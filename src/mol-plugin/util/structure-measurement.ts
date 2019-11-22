/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '../../mol-model/structure';
import { PluginContext } from '../context';
import { StateSelection, StateTransform } from '../../mol-state';
import { StateTransforms } from '../state/transforms';
import { PluginCommands } from '../command';

export { StructureMeasurementManager }

const MeasurementGroupTag = 'measurement-group';

class StructureMeasurementManager {
    private getGroup() {
        const state = this.context.state.dataState;
        const groupRef = StateSelection.findTagInSubtree(state.tree, StateTransform.RootRef, MeasurementGroupTag);
        const builder = this.context.state.dataState.build();

        if (groupRef) return builder.to(groupRef);
        return builder.toRoot().group(StateTransforms.Misc.CreateGroup, { label: `Measurements` }, { tags: MeasurementGroupTag });
    }

    async addDistance(a: StructureElement.Loci, b: StructureElement.Loci) {
        const cellA = this.context.helpers.substructureParent.get(a.structure);
        const cellB = this.context.helpers.substructureParent.get(b.structure);

        if (!cellA || !cellB) return;

        const update = this.getGroup();

        console.log({ cellA, cellB });

        const state = this.context.state.dataState;
        await PluginCommands.State.Update.dispatch(this.context, { state, tree: update, options: { doNotLogTiming: true } });
    }

    constructor(private context: PluginContext) {

    }
}