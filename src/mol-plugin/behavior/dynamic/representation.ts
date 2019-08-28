/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MarkerAction } from '../../../mol-util/marker-action';
import { PluginContext } from '../../../mol-plugin/context';
import { PluginStateObject as SO } from '../../state/objects';
import { lociLabel } from '../../../mol-theme/label';
import { PluginBehavior } from '../behavior';
import { Interactivity } from '../../util/interactivity';
import { StateTreeSpine } from '../../../mol-state/tree/spine';
import { StateSelection } from '../../../mol-state';

export const HighlightLoci = PluginBehavior.create({
    name: 'representation-highlight-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler {
        private lociMarkProvider = (interactionLoci: Interactivity.Loci, action: MarkerAction) => {
            if (!this.ctx.canvas3d) return;
            this.ctx.canvas3d.mark({ loci: interactionLoci.loci }, action)
        }
        register() {
            this.ctx.interactivity.lociHighlights.addProvider(this.lociMarkProvider)
        }
        unregister() {
            this.ctx.interactivity.lociHighlights.removeProvider(this.lociMarkProvider)
        }
    },
    display: { name: 'Highlight Loci on Canvas' }
});

export const SelectLoci = PluginBehavior.create({
    name: 'representation-select-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler {
        private spine: StateTreeSpine.Impl
        private lociMarkProvider = (interactionLoci: Interactivity.Loci, action: MarkerAction) => {
            if (!this.ctx.canvas3d) return;
            this.ctx.canvas3d.mark({ loci: interactionLoci.loci }, action)
        }
        private applySelectMark(ref: string) {
            const cell = this.ctx.state.dataState.cells.get(ref)
            if (cell && SO.isRepresentation3D(cell.obj)) {
                this.spine.current = cell
                const so = this.spine.getRootOfType(SO.Molecule.Structure)
                if (so) {
                    const loci = this.ctx.helpers.structureSelectionManager.get(so.data)
                    this.lociMarkProvider({ loci }, MarkerAction.Select)
                }
            }
        }
        register() {
            this.ctx.interactivity.lociSelections.addProvider(this.lociMarkProvider)

            this.subscribeObservable(this.ctx.events.state.object.created, ({ ref }) => this.applySelectMark(ref));

            // re-apply select-mark to all representation of an updated structure
            this.subscribeObservable(this.ctx.events.state.object.updated, ({ ref }) => {
                const cell = this.ctx.state.dataState.cells.get(ref)
                if (cell && SO.Molecule.Structure.is(cell.obj)) {
                    const reprs = this.ctx.state.dataState.select(StateSelection.Generators.ofType(SO.Molecule.Structure.Representation3D, ref))
                    for (const repr of reprs) this.applySelectMark(repr.transform.ref)
                }
            });
        }
        unregister() {
            this.ctx.interactivity.lociSelections.removeProvider(this.lociMarkProvider)
        }
        constructor(ctx: PluginContext, params: {}) {
            super(ctx, params)
            this.spine = new StateTreeSpine.Impl(ctx.state.dataState.cells)
        }
    },
    display: { name: 'Select Loci on Canvas' }
});

export const DefaultLociLabelProvider = PluginBehavior.create({
    name: 'default-loci-label-provider',
    category: 'interaction',
    ctor: class implements PluginBehavior<undefined> {
        private f = lociLabel;
        register() { this.ctx.lociLabels.addProvider(this.f); }
        unregister() { this.ctx.lociLabels.removeProvider(this.f); }
        constructor(protected ctx: PluginContext) { }
    },
    display: { name: 'Provide Default Loci Label' }
});
