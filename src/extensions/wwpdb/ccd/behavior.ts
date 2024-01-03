/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { PluginBehavior } from '../../../mol-plugin/behavior/behavior';
import { ChemicalComponentPreset, ChemicalCompontentTrajectoryHierarchyPreset } from './representation';

export const wwPDBChemicalComponentDictionary = PluginBehavior.create<{ }>({
    name: 'wwpdb-chemical-component-dictionary',
    category: 'representation',
    display: {
        name: 'wwPDB Chemical Compontent Dictionary',
        description: 'Custom representation for data loaded from the CCD.'
    },
    ctor: class extends PluginBehavior.Handler<{ }> {
        register(): void {
            this.ctx.builders.structure.hierarchy.registerPreset(ChemicalCompontentTrajectoryHierarchyPreset);
            this.ctx.builders.structure.representation.registerPreset(ChemicalComponentPreset);
        }

        update() {
            return false;
        }

        unregister() {
            this.ctx.builders.structure.hierarchy.unregisterPreset(ChemicalCompontentTrajectoryHierarchyPreset);
            this.ctx.builders.structure.representation.unregisterPreset(ChemicalComponentPreset);
        }
    },
    params: () => ({ })
});