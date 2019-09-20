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
import { ButtonsType, ModifiersKeys } from '../../../mol-util/input/input-observer';
import { Binding } from '../../../mol-util/binding';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

const B = ButtonsType
const M = ModifiersKeys
const Trigger = Binding.Trigger

//

const DefaultHighlightLociBindings = {
    hoverHighlightOnly: Binding(Trigger(B.Flag.None), 'Highlight hovered element using ${trigger}'),
    hoverHighlightOnlyExtend: Binding(Trigger(B.Flag.None, M.create({ shift: true })), 'Extend highlight from selected to hovered element along polymer using ${trigger}'),
}
const HighlightLociParams = {
    bindings: PD.Value(DefaultHighlightLociBindings, { isHidden: true }),
}
type HighlightLociProps = PD.Values<typeof HighlightLociParams>

export const HighlightLoci = PluginBehavior.create({
    name: 'representation-highlight-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<HighlightLociProps> {
        private lociMarkProvider = (interactionLoci: Interactivity.Loci, action: MarkerAction) => {
            if (!this.ctx.canvas3d) return;
            this.ctx.canvas3d.mark({ loci: interactionLoci.loci }, action)
        }
        register() {
            this.subscribeObservable(this.ctx.behaviors.interaction.hover, ({ current, buttons, modifiers }) => {
                if (!this.ctx.canvas3d) return

                if (Binding.match(this.params.bindings.hoverHighlightOnly, buttons, modifiers)) {
                    this.ctx.interactivity.lociHighlights.highlightOnly(current)
                }

                if (Binding.match(this.params.bindings.hoverHighlightOnlyExtend, buttons, modifiers)) {
                    this.ctx.interactivity.lociHighlights.highlightOnlyExtend(current)
                }
            });
            this.ctx.interactivity.lociHighlights.addProvider(this.lociMarkProvider)
        }
        unregister() {
            this.ctx.interactivity.lociHighlights.removeProvider(this.lociMarkProvider)
        }
    },
    params: () => HighlightLociParams,
    display: { name: 'Highlight Loci on Canvas' }
});

//

const DefaultSelectLociBindings = {
    clickSelect: Binding.Empty,
    clickSelectExtend: Binding(Trigger(B.Flag.Primary, M.create({ shift: true })), 'Extend selection to clicked element along polymer using ${trigger}.'),
    clickSelectOnly: Binding(Trigger(B.Flag.Secondary, M.create({ control: true })), 'Select only the clicked element using ${trigger}.'),
    clickSelectToggle: Binding(Trigger(B.Flag.Primary, M.create({ control: true })), 'Toggle selection of clicked element using ${trigger}.'),
    clickDeselect: Binding.Empty,
}
const SelectLociParams = {
    bindings: PD.Value(DefaultSelectLociBindings, { isHidden: true }),
}
type SelectLociProps = PD.Values<typeof SelectLociParams>

export const SelectLoci = PluginBehavior.create({
    name: 'representation-select-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<SelectLociProps> {
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
            this.subscribeObservable(this.ctx.behaviors.interaction.click, ({ current, buttons, modifiers }) => {
                if (!this.ctx.canvas3d) return

                if (Binding.match(this.params.bindings.clickSelect, buttons, modifiers)) {
                    this.ctx.interactivity.lociSelects.select(current)
                }

                if (Binding.match(this.params.bindings.clickSelectExtend, buttons, modifiers)) {
                    this.ctx.interactivity.lociSelects.selectExtend(current)
                }

                if (Binding.match(this.params.bindings.clickSelectOnly, buttons, modifiers)) {
                    this.ctx.interactivity.lociSelects.selectOnly(current)
                }

                if (Binding.match(this.params.bindings.clickSelectToggle, buttons, modifiers)) {
                    this.ctx.interactivity.lociSelects.selectToggle(current)
                }

                if (Binding.match(this.params.bindings.clickDeselect, buttons, modifiers)) {
                    this.ctx.interactivity.lociSelects.deselect(current)
                }
            });
            this.ctx.interactivity.lociSelects.addProvider(this.lociMarkProvider)

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
            this.ctx.interactivity.lociSelects.removeProvider(this.lociMarkProvider)
        }
        constructor(ctx: PluginContext, params: SelectLociProps) {
            super(ctx, params)
            this.spine = new StateTreeSpine.Impl(ctx.state.dataState.cells)
        }
    },
    params: () => SelectLociParams,
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
