/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MarkerAction } from '../../../mol-util/marker-action';
import { PluginContext } from '../../../mol-plugin/context';
import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { lociLabel } from '../../../mol-theme/label';
import { PluginBehavior } from '../behavior';
import { StateTreeSpine } from '../../../mol-state/tree/spine';
import { StateSelection } from '../../../mol-state';
import { ButtonsType, ModifiersKeys } from '../../../mol-util/input/input-observer';
import { Binding } from '../../../mol-util/binding';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Structure, StructureElement, StructureProperties } from '../../../mol-model/structure';
import { arrayMax } from '../../../mol-util/array';
import { Representation } from '../../../mol-repr/representation';
import { LociLabel } from '../../../mol-plugin-state/manager/loci-label';

const B = ButtonsType;
const M = ModifiersKeys;
const Trigger = Binding.Trigger;

//

const DefaultHighlightLociBindings = {
    hoverHighlightOnly: Binding([Trigger(B.Flag.None)], 'Highlight', 'Hover element using ${triggers}'),
    hoverHighlightOnlyExtend: Binding([Trigger(B.Flag.None, M.create({ shift: true }))], 'Extend highlight', 'From selected to hovered element along polymer using ${triggers}'),
};
const HighlightLociParams = {
    bindings: PD.Value(DefaultHighlightLociBindings, { isHidden: true }),
};
type HighlightLociProps = PD.Values<typeof HighlightLociParams>

export const HighlightLoci = PluginBehavior.create({
    name: 'representation-highlight-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<HighlightLociProps> {
        private lociMarkProvider = (interactionLoci: Representation.Loci, action: MarkerAction) => {
            if (!this.ctx.canvas3d) return;
            this.ctx.canvas3d.mark({ loci: interactionLoci.loci }, action);
        }
        register() {
            this.subscribeObservable(this.ctx.behaviors.interaction.hover, ({ current, buttons, modifiers }) => {
                if (!this.ctx.canvas3d || this.ctx.isBusy) return;
                let matched = false;

                if (Binding.match(this.params.bindings.hoverHighlightOnly, buttons, modifiers)) {
                    this.ctx.managers.interactivity.lociHighlights.highlightOnly(current);
                    matched = true;
                }

                if (Binding.match(this.params.bindings.hoverHighlightOnlyExtend, buttons, modifiers)) {
                    this.ctx.managers.interactivity.lociHighlights.highlightOnlyExtend(current);
                    matched = true;
                }

                if (!matched) {
                    this.ctx.managers.interactivity.lociHighlights.highlightOnly({ repr: current.repr, loci: EmptyLoci });
                }
            });
            this.ctx.managers.interactivity.lociHighlights.addProvider(this.lociMarkProvider);
        }
        unregister() {
            this.ctx.managers.interactivity.lociHighlights.removeProvider(this.lociMarkProvider);
        }
    },
    params: () => HighlightLociParams,
    display: { name: 'Highlight Loci on Canvas' }
});

//

const DefaultSelectLociBindings = {
    clickSelect: Binding.Empty,
    clickToggleExtend: Binding([Trigger(B.Flag.Primary, M.create({ shift: true }))], 'Toggle extended selection', '${triggers} to extend selection along polymer'),
    clickSelectOnly: Binding.Empty,
    clickToggle: Binding([Trigger(B.Flag.Primary, M.create())], 'Toggle selection', '${triggers} on element'),
    clickDeselect: Binding.Empty,
    clickDeselectAllOnEmpty: Binding([Trigger(B.Flag.Primary, M.create())], 'Deselect all', 'Click on nothing using ${triggers}'),
};
const SelectLociParams = {
    bindings: PD.Value(DefaultSelectLociBindings, { isHidden: true }),
};
type SelectLociProps = PD.Values<typeof SelectLociParams>

export const SelectLoci = PluginBehavior.create({
    name: 'representation-select-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<SelectLociProps> {
        private spine: StateTreeSpine.Impl
        private lociMarkProvider = (reprLoci: Representation.Loci, action: MarkerAction) => {
            if (!this.ctx.canvas3d) return;
            this.ctx.canvas3d.mark({ loci: reprLoci.loci }, action);
        }
        private applySelectMark(ref: string, clear?: boolean) {
            const cell = this.ctx.state.data.cells.get(ref);
            if (cell && SO.isRepresentation3D(cell.obj)) {
                this.spine.current = cell;
                const so = this.spine.getRootOfType(SO.Molecule.Structure);
                if (so) {
                    if (clear) {
                        this.lociMarkProvider({ loci: Structure.Loci(so.data) }, MarkerAction.Deselect);
                    }
                    const loci = this.ctx.managers.structure.selection.getLoci(so.data);
                    this.lociMarkProvider({ loci }, MarkerAction.Select);
                }
            }
        }
        register() {
            const lociIsEmpty = (current: Representation.Loci) => Loci.isEmpty(current.loci);
            const lociIsNotEmpty = (current: Representation.Loci) => !Loci.isEmpty(current.loci);

            const actions: [keyof typeof DefaultSelectLociBindings, (current: Representation.Loci) => void, ((current: Representation.Loci) => boolean) | undefined][] = [
                ['clickSelect', current => this.ctx.managers.interactivity.lociSelects.select(current), lociIsNotEmpty],
                ['clickToggle', current => this.ctx.managers.interactivity.lociSelects.toggle(current), lociIsNotEmpty],
                ['clickToggleExtend', current => this.ctx.managers.interactivity.lociSelects.toggleExtend(current), lociIsNotEmpty],
                ['clickSelectOnly', current => this.ctx.managers.interactivity.lociSelects.selectOnly(current), lociIsNotEmpty],
                ['clickDeselect', current => this.ctx.managers.interactivity.lociSelects.deselect(current), lociIsNotEmpty],
                ['clickDeselectAllOnEmpty', () => this.ctx.managers.interactivity.lociSelects.deselectAll(), lociIsEmpty],
            ];

            // sort the action so that the ones with more modifiers trigger sooner.
            actions.sort((a, b) => {
                const x = this.params.bindings[a[0]], y = this.params.bindings[b[0]];
                const k = x.triggers.length === 0 ? 0 : arrayMax(x.triggers.map(t => M.size(t.modifiers)));
                const l = y.triggers.length === 0 ? 0 : arrayMax(y.triggers.map(t => M.size(t.modifiers)));
                return l - k;
            });

            this.subscribeObservable(this.ctx.behaviors.interaction.click, ({ current, button, modifiers }) => {
                if (!this.ctx.canvas3d || this.ctx.isBusy || !this.ctx.selectionMode) return;

                // only trigger the 1st action that matches
                for (const [binding, action, condition] of actions) {
                    if (Binding.match(this.params.bindings[binding], button, modifiers) && (!condition || condition(current))) {
                        action(current);
                        break;
                    }
                }
            });
            this.ctx.managers.interactivity.lociSelects.addProvider(this.lociMarkProvider);

            this.subscribeObservable(this.ctx.state.events.object.created, ({ ref }) => this.applySelectMark(ref));

            // re-apply select-mark to all representation of an updated structure
            this.subscribeObservable(this.ctx.state.events.object.updated, ({ ref }) => {
                const cell = this.ctx.state.data.cells.get(ref);
                if (cell && SO.Molecule.Structure.is(cell.obj)) {
                    const reprs = this.ctx.state.data.select(StateSelection.Generators.ofType(SO.Molecule.Structure.Representation3D, ref));
                    for (const repr of reprs) this.applySelectMark(repr.transform.ref, true);
                }
            });
        }
        unregister() {
            this.ctx.managers.interactivity.lociSelects.removeProvider(this.lociMarkProvider);
        }
        constructor(ctx: PluginContext, params: SelectLociProps) {
            super(ctx, params);
            this.spine = new StateTreeSpine.Impl(ctx.state.data.cells);
        }
    },
    params: () => SelectLociParams,
    display: { name: 'Select Loci on Canvas' }
});

//

export const DefaultLociLabelProvider = PluginBehavior.create({
    name: 'default-loci-label-provider',
    category: 'interaction',
    ctor: class implements PluginBehavior<undefined> {
        private f = {
            label: (loci: Loci) => {
                const label: string[] = [];
                if (StructureElement.Loci.is(loci) && loci.elements.length === 1) {
                    const { unit: u } = loci.elements[0];
                    const l = StructureElement.Location.create(loci.structure, u, u.elements[0]);
                    const name = StructureProperties.entity.pdbx_description(l).join(', ');
                    label.push(name);
                }
                label.push(lociLabel(loci));
                return label.filter(l => !!l).join('</br>');
            },
            group: (label: LociLabel) => label.toString().replace(/Model [0-9]+/g, 'Models'),
            priority: 100
        };
        register() { this.ctx.managers.lociLabels.addProvider(this.f); }
        unregister() { this.ctx.managers.lociLabels.removeProvider(this.f); }
        constructor(protected ctx: PluginContext) { }
    },
    display: { name: 'Provide Default Loci Label' }
});

//

const DefaultFocusLociBindings = {
    clickFocus: Binding([
        Trigger(B.Flag.Primary, M.create()),
    ], 'Representation Focus', 'Click element using ${triggers}'),
    clickFocusAdd: Binding([
        Trigger(B.Flag.Primary, M.create({ shift: true })),
    ], 'Representation Focus Add', 'Click element using ${triggers}'),
    clickFocusSelectMode: Binding([
        // default is empty
    ], 'Representation Focus', 'Click element using ${triggers}'),
    clickFocusAddSelectMode: Binding([
        // default is empty
    ], 'Representation Focus Add', 'Click element using ${triggers}'),
};
const FocusLociParams = {
    bindings: PD.Value(DefaultFocusLociBindings, { isHidden: true }),
};
type FocusLociProps = PD.Values<typeof FocusLociParams>

export const FocusLoci = PluginBehavior.create<FocusLociProps>({
    name: 'representation-focus-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<FocusLociProps> {
        register(): void {
            this.subscribeObservable(this.ctx.behaviors.interaction.click, ({ current, button, modifiers }) => {
                const { clickFocus, clickFocusAdd, clickFocusSelectMode, clickFocusAddSelectMode } = this.params.bindings;

                // only apply structure focus for appropriate granularity
                const { granularity } = this.ctx.managers.interactivity.props;
                if (granularity !== 'residue' && granularity !== 'element') return;

                const binding = this.ctx.selectionMode ? clickFocusSelectMode : clickFocus;
                const matched = Binding.match(binding, button, modifiers);
                const bindingAdd = this.ctx.selectionMode ? clickFocusAddSelectMode : clickFocusAdd;
                const matchedAdd = Binding.match(bindingAdd, button, modifiers);
                if (!matched && !matchedAdd) return;

                const loci = Loci.normalize(current.loci, 'residue');
                const entry = this.ctx.managers.structure.focus.current;
                if (entry && Loci.areEqual(entry.loci, loci)) {
                    this.ctx.managers.structure.focus.clear();
                } else {
                    if (matched) {
                        this.ctx.managers.structure.focus.setFromLoci(loci);
                    } else {
                        this.ctx.managers.structure.focus.addFromLoci(loci);

                        // focus-add is not handled in camera behavior, doing it here
                        const current = this.ctx.managers.structure.focus.current?.loci;
                        if (current) this.ctx.managers.camera.focusLoci(current);
                    }
                }
            });
        }
    },
    params: () => FocusLociParams,
    display: { name: 'Representation Focus Loci on Canvas' }
});