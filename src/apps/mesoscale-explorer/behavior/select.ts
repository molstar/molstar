/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { EveryLoci, Loci } from '../../../mol-model/loci';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginBehavior } from '../../../mol-plugin/behavior';
import { ButtonsType, ModifiersKeys } from '../../../mol-util/input/input-observer';
import { Binding } from '../../../mol-util/binding';
import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { Structure, StructureElement } from '../../../mol-model/structure';
import { StateSelection } from '../../../mol-state';
import { StateTreeSpine } from '../../../mol-state/tree/spine';
import { Representation } from '../../../mol-repr/representation';
import { MarkerAction } from '../../../mol-util/marker-action';
import { PluginContext } from '../../../mol-plugin/context';

const B = ButtonsType;
const M = ModifiersKeys;
const Trigger = Binding.Trigger;

const DefaultMesoSelectLociBindings = {
    clickToggleSelect: Binding([
        Trigger(B.Flag.Primary, M.create({ shift: true })),
        Trigger(B.Flag.Primary, M.create({ control: true })),
    ], 'Toggle select', 'Click element using ${triggers}'),
    hoverHighlightOnly: Binding([
        Trigger(B.Flag.None, M.create({ shift: true })),
        Trigger(B.Flag.None, M.create({ control: true })),
    ], 'Highlight', 'Hover element using ${triggers}'),
};
const MesoSelectLociParams = {
    bindings: PD.Value(DefaultMesoSelectLociBindings, { isHidden: true }),
};
type MesoSelectLociProps = PD.Values<typeof MesoSelectLociParams>

export const MesoSelectLoci = PluginBehavior.create<MesoSelectLociProps>({
    name: 'camera-meso-select-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<MesoSelectLociProps> {
        private spine: StateTreeSpine.Impl;
        private lociMarkProvider = (interactionLoci: Representation.Loci, action: MarkerAction) => {
            if (!this.ctx.canvas3d) return;
            this.ctx.canvas3d.mark(interactionLoci, action);
        };
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
        register(): void {
            this.subscribeObservable(this.ctx.behaviors.interaction.click, ({ current, button, modifiers }) => {
                if (!this.ctx.canvas3d || this.ctx.isBusy) return;

                const { clickToggleSelect } = this.params.bindings;
                if (Binding.match(clickToggleSelect, button, modifiers)) {
                    if (Loci.isEmpty(current.loci)) {
                        this.ctx.managers.interactivity.lociSelects.deselectAll();
                        return;
                    }

                    const loci = Loci.normalize(current.loci, modifiers.control ? 'entity' : 'chain');
                    this.ctx.managers.interactivity.lociSelects.toggle({ loci }, false);
                }
            });
            this.ctx.managers.interactivity.lociSelects.addProvider(this.lociMarkProvider);

            this.subscribeObservable(this.ctx.behaviors.interaction.hover, ({ current, button, modifiers }) => {
                if (!this.ctx.canvas3d || this.ctx.isBusy) return;

                const pointerLock = !!this.ctx.canvas3dContext?.input.pointerLock;
                const { hoverHighlightOnly } = this.params.bindings;

                if (!pointerLock && Binding.match(hoverHighlightOnly, button, modifiers)) {
                    if (Loci.isEmpty(current.loci)) {
                        this.ctx.managers.interactivity.lociHighlights.clearHighlights();
                        return;
                    }

                    if (modifiers.control) {
                        this.ctx.managers.interactivity.lociHighlights.highlightOnly({ repr: current.repr, loci: EveryLoci }, false);
                    } else {
                        const loci = Loci.normalize(current.loci, 'chain');
                        this.ctx.managers.interactivity.lociHighlights.highlightOnly({ repr: current.repr, loci }, false);
                    }
                }

                if (Loci.isEmpty(current.loci)) {
                    this.ctx.behaviors.labels.highlight.next({ labels: [] });
                } else {
                    const labels: string[] = [];
                    if (StructureElement.Loci.is(current.loci)) {
                        const cell = this.ctx.helpers.substructureParent.get(current.loci.structure);
                        labels.push(cell?.obj?.label || 'Unknown');
                    }
                    this.ctx.behaviors.labels.highlight.next({ labels });
                }
            });
            this.ctx.managers.interactivity.lociHighlights.addProvider(this.lociMarkProvider);

            let dimDisabled = false;

            this.subscribeObservable(this.ctx.behaviors.interaction.keyReleased, ({ code, modifiers }) => {
                if (!this.ctx.canvas3d) return;

                if ((code.startsWith('Shift') && !modifiers.control) || (code.startsWith('Control') && !modifiers.shift)) {
                    if (dimDisabled) {
                        dimDisabled = false;
                        this.ctx.canvas3d?.setProps({ renderer: { dimStrength: 1 } }, true);
                    }
                    this.ctx.managers.interactivity.lociHighlights.clearHighlights();
                }
            });

            this.subscribeObservable(this.ctx.behaviors.interaction.key, ({ modifiers }) => {
                if (!this.ctx.canvas3d) return;

                if (!dimDisabled && modifiers.control && modifiers.shift) {
                    dimDisabled = true;
                    this.ctx.canvas3d?.setProps({ renderer: { dimStrength: 0 } });
                }
            });

            this.subscribeObservable(this.ctx.state.events.object.created, ({ ref }) => this.applySelectMark(ref));

            // re-apply select-mark to all representation of an updated structure
            this.subscribeObservable(this.ctx.state.events.object.updated, ({ ref, obj, oldObj, oldData, action }) => {
                const cell = this.ctx.state.data.cells.get(ref);
                if (cell && SO.Molecule.Structure.is(cell.obj)) {
                    const structure: Structure = obj.data;
                    const oldStructure: Structure | undefined = action === 'recreate' ? oldObj?.data :
                        action === 'in-place' ? oldData : undefined;
                    if (oldStructure &&
                        Structure.areEquivalent(structure, oldStructure) &&
                        Structure.areHierarchiesEqual(structure, oldStructure)) return;

                    const reprs = this.ctx.state.data.select(StateSelection.children(ref).ofType(SO.Molecule.Structure.Representation3D));
                    for (const repr of reprs) this.applySelectMark(repr.transform.ref, true);
                }
            });
        }
        unregister() {
            this.ctx.managers.interactivity.lociSelects.removeProvider(this.lociMarkProvider);
            this.ctx.managers.interactivity.lociHighlights.removeProvider(this.lociMarkProvider);
        }
        constructor(ctx: PluginContext, params: MesoSelectLociProps) {
            super(ctx, params);
            this.spine = new StateTreeSpine.Impl(ctx.state.data.cells);
        }
    },
    params: () => MesoSelectLociParams,
    display: { name: 'Camera Meso Select Loci on Canvas' }
});
