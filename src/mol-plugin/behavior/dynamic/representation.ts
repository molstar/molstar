/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MarkerAction } from '../../../mol-geo/geometry/marker-data';
import { EmptyLoci } from '../../../mol-model/loci';
import { StructureElement } from '../../../mol-model/structure';
import { PluginContext } from '../../../mol-plugin/context';
import { Representation } from '../../../mol-repr/representation';
import { labelFirst } from '../../../mol-theme/label';
import { ButtonsType } from '../../../mol-util/input/input-observer';
import { PluginBehavior } from '../behavior';

export const HighlightLoci = PluginBehavior.create({
    name: 'representation-highlight-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler {
        register(): void {
            let prev: Representation.Loci = { loci: EmptyLoci, repr: void 0 };
            const sel = this.ctx.helpers.structureSelection;

            this.subscribeObservable(this.ctx.behaviors.canvas3d.highlight, ({ current, modifiers }) => {
                if (!this.ctx.canvas3d) return;

                if (StructureElement.isLoci(current.loci)) {
                    let loci: StructureElement.Loci = current.loci;
                    if (modifiers && modifiers.shift) {
                        loci = sel.tryGetRange(loci) || loci;
                    }

                    this.ctx.canvas3d.mark(prev, MarkerAction.RemoveHighlight);
                    const toHighlight = { loci, repr: current.repr };
                    this.ctx.canvas3d.mark(toHighlight, MarkerAction.Highlight);
                    prev = toHighlight;
                } else {
                    if (!Representation.Loci.areEqual(prev, current)) {
                        this.ctx.canvas3d.mark(prev, MarkerAction.RemoveHighlight);
                        this.ctx.canvas3d.mark(current, MarkerAction.Highlight);
                        prev = current;
                    }
                }

            });
        }
    },
    display: { name: 'Highlight Loci on Canvas' }
});

export const SelectLoci = PluginBehavior.create({
    name: 'representation-select-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler {
        register(): void {
            const sel = this.ctx.helpers.structureSelection;

            const toggleSel = (current: Representation.Loci<StructureElement.Loci>) => {
                if (sel.has(current.loci)) {
                    sel.remove(current.loci);
                    this.ctx.canvas3d.mark(current, MarkerAction.Deselect);
                } else {
                    sel.add(current.loci);
                    this.ctx.canvas3d.mark(current, MarkerAction.Select);
                }
            }

            this.subscribeObservable(this.ctx.behaviors.canvas3d.click, ({ current, buttons, modifiers }) => {
                if (!this.ctx.canvas3d) return;

                if (current.loci.kind === 'empty-loci') {
                    if (modifiers.control && buttons === ButtonsType.Flag.Secondary) {
                        // clear the selection on Ctrl + Right-Click on empty
                        const sels = sel.clear();
                        for (const s of sels) this.ctx.canvas3d.mark({ loci: s }, MarkerAction.Deselect);
                    }
                } else if (StructureElement.isLoci(current.loci)) {
                    if (modifiers.control && buttons === ButtonsType.Flag.Secondary) {
                        // select only the current element on Ctrl + Right-Click
                        const old = sel.get(current.loci.structure);
                        this.ctx.canvas3d.mark({ loci: old }, MarkerAction.Deselect);
                        sel.set(current.loci);
                        this.ctx.canvas3d.mark(current, MarkerAction.Select);
                    } else if (modifiers.control && buttons === ButtonsType.Flag.Primary) {
                        // toggle current element on Ctrl + Left-Click
                        toggleSel(current as Representation.Loci<StructureElement.Loci>);
                    } else if (modifiers.shift && buttons === ButtonsType.Flag.Primary) {
                        // try to extend sequence on Shift + Left-Click
                        let loci: StructureElement.Loci = current.loci;
                        if (modifiers && modifiers.shift) {
                            loci = sel.tryGetRange(loci) || loci;
                        }
                        toggleSel({ loci, repr: current.repr });
                    }
                } else {
                    if (!ButtonsType.has(buttons, ButtonsType.Flag.Secondary)) return;
                    this.ctx.canvas3d.mark(current, MarkerAction.Toggle);
                }
            });
        }
    },
    display: { name: 'Select Loci on Canvas' }
});

export const DefaultLociLabelProvider = PluginBehavior.create({
    name: 'default-loci-label-provider',
    category: 'interaction',
    ctor: class implements PluginBehavior<undefined> {
        private f = labelFirst;
        register(): void { this.ctx.lociLabels.addProvider(this.f); }
        unregister() { this.ctx.lociLabels.removeProvider(this.f); }
        constructor(protected ctx: PluginContext) { }
    },
    display: { name: 'Provide Default Loci Label' }
});
