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
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Structure } from '../../../mol-model/structure';
import { arrayMax } from '../../../mol-util/array';

const B = ButtonsType
const M = ModifiersKeys
const Trigger = Binding.Trigger

//

const DefaultHighlightLociBindings = {
    hoverHighlightOnly: Binding([Trigger(B.Flag.None)], 'Highlight hovered element using ${triggers}'),
    hoverHighlightOnlyExtend: Binding([Trigger(B.Flag.None, M.create({ shift: true }))], 'Extend highlight from selected to hovered element along polymer using ${triggers}'),
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
            this.ctx.canvas3d.mark({ loci: interactionLoci.loci, repr: interactionLoci.passRepresentation ? interactionLoci.repr : void 0 }, action)
        }
        register() {
            this.subscribeObservable(this.ctx.behaviors.interaction.hover, ({ current, buttons, modifiers }) => {
                if (!this.ctx.canvas3d) return
                let matched = false

                if (Binding.match(this.params.bindings.hoverHighlightOnly, buttons, modifiers)) {
                    this.ctx.interactivity.lociHighlights.highlightOnly(current)
                    matched = true
                }

                if (Binding.match(this.params.bindings.hoverHighlightOnlyExtend, buttons, modifiers)) {
                    this.ctx.interactivity.lociHighlights.highlightOnlyExtend(current)
                    matched = true
                }

                if (!matched) {
                    this.ctx.interactivity.lociHighlights.highlightOnly({ repr: current.repr, loci: EmptyLoci })
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
    clickToggleExtend: Binding([Trigger(B.Flag.Primary, M.create({ shift: true }))], 'Toggle selection of last click element extended along polymer using ${triggers}.'),
    clickSelectOnly: Binding.Empty,
    clickToggle: Binding([Trigger(B.Flag.Primary, M.create())], 'Toggle selection of clicked element using ${triggers}.'),
    clickDeselect: Binding.Empty,
    clickDeselectAllOnEmpty: Binding([Trigger(B.Flag.Primary, M.create())], 'Deselect all when clicking on nothing using ${triggers}.'),
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
            this.ctx.canvas3d.mark({ loci: interactionLoci.loci, repr: interactionLoci.passRepresentation ? interactionLoci.repr : void 0 }, action)
        }
        private applySelectMark(ref: string, clear?: boolean) {
            const cell = this.ctx.state.dataState.cells.get(ref)
            if (cell && SO.isRepresentation3D(cell.obj)) {
                this.spine.current = cell
                const so = this.spine.getRootOfType(SO.Molecule.Structure)
                if (so) {
                    if (clear) {
                        this.lociMarkProvider({ loci: Structure.Loci(so.data) }, MarkerAction.Deselect)
                    }
                    const loci = this.ctx.helpers.structureSelectionManager.get(so.data)
                    this.lociMarkProvider({ loci }, MarkerAction.Select)
                }
            }
        }
        register() {
            const lociIsEmpty = (current: Interactivity.Loci) => Loci.isEmpty(current.loci)
            const lociIsNotEmpty = (current: Interactivity.Loci) => !Loci.isEmpty(current.loci)

            const actions: [keyof typeof DefaultSelectLociBindings, (current: Interactivity.Loci) => void, ((current: Interactivity.Loci) => boolean) | undefined][] = [
                ['clickSelect', current => this.ctx.interactivity.lociSelects.select(current), lociIsNotEmpty],
                ['clickToggle', current => this.ctx.interactivity.lociSelects.toggle(current), lociIsNotEmpty],
                ['clickToggleExtend', current => this.ctx.interactivity.lociSelects.toggleExtend(current), lociIsNotEmpty],
                ['clickSelectOnly', current => this.ctx.interactivity.lociSelects.selectOnly(current), lociIsNotEmpty],
                ['clickDeselect', current => this.ctx.interactivity.lociSelects.deselect(current), lociIsNotEmpty],
                ['clickDeselectAllOnEmpty', () => this.ctx.interactivity.lociSelects.deselectAll(), lociIsEmpty],
            ];

            // sort the action so that the ones with more modifiers trigger sooner.
            actions.sort((a, b) => {
                const x = this.params.bindings[a[0]], y = this.params.bindings[b[0]];
                const k = x.triggers.length === 0 ? 0 : arrayMax(x.triggers.map(t => M.size(t.modifiers)));
                const l = y.triggers.length === 0 ? 0 : arrayMax(y.triggers.map(t => M.size(t.modifiers)));
                return l - k;
            })

            this.subscribeObservable(this.ctx.behaviors.interaction.click, ({ current, button, modifiers }) => {
                if (!this.ctx.canvas3d) return

                // only trigger the 1st action that matches
                for (const [binding, action, condition] of actions) {
                    if (Binding.match(this.params.bindings[binding], button, modifiers) && (!condition || condition(current))) {
                        action(current);
                        break;
                    }
                }
            });
            this.ctx.interactivity.lociSelects.addProvider(this.lociMarkProvider)

            this.subscribeObservable(this.ctx.events.state.object.created, ({ ref }) => this.applySelectMark(ref));

            // re-apply select-mark to all representation of an updated structure
            this.subscribeObservable(this.ctx.events.state.object.updated, ({ ref }) => {
                const cell = this.ctx.state.dataState.cells.get(ref)
                if (cell && SO.Molecule.Structure.is(cell.obj)) {
                    const reprs = this.ctx.state.dataState.select(StateSelection.Generators.ofType(SO.Molecule.Structure.Representation3D, ref))
                    for (const repr of reprs) this.applySelectMark(repr.transform.ref, true)
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
        private f = (loci: Loci) => lociLabel(loci);
        register() { this.ctx.lociLabels.addProvider(this.f); }
        unregister() { this.ctx.lociLabels.removeProvider(this.f); }
        constructor(protected ctx: PluginContext) { }
    },
    display: { name: 'Provide Default Loci Label' }
});
