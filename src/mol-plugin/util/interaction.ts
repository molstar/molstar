/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Loci as ModelLoci, EmptyLoci } from '../../mol-model/loci';
import { ModifiersKeys, ButtonsType } from '../../mol-util/input/input-observer';
import { Representation } from '../../mol-repr/representation';
import { StructureElement } from '../../mol-model/structure';
import { MarkerAction } from '../../mol-util/marker-action';
import { StructureElementSelectionManager } from './structure-element-selection';
import { PluginContext } from '../context';

export namespace Interaction {
    export interface Loci<T extends ModelLoci = ModelLoci> { loci: T, repr?: Representation.Any }

    export namespace Loci {
        export function areEqual(a: Loci, b: Loci) {
            return a.repr === b.repr && ModelLoci.areEqual(a.loci, b.loci);
        }
        export const Empty: Loci = { loci: EmptyLoci };
    }

    export interface HighlightEvent { current: Loci, modifiers?: ModifiersKeys }
    export interface ClickEvent { current: Loci, buttons: ButtonsType, modifiers: ModifiersKeys }

    export type LociMarkProvider = (loci: Loci, action: MarkerAction) => void

    export abstract class LociMarkManager<MarkEvent extends any> {
        protected providers: LociMarkProvider[] = [];
        protected sel: StructureElementSelectionManager

        addProvider(provider: LociMarkProvider) {
            this.providers.push(provider);
        }

        removeProvider(provider: LociMarkProvider) {
            this.providers = this.providers.filter(p => p !== provider);
            // TODO clear, then re-apply remaining providers
        }

        toggleSel(current: Loci<ModelLoci>) {
            if (this.sel.has(current.loci)) {
                this.sel.remove(current.loci);
                this.mark(current, MarkerAction.Deselect);
            } else {
                this.sel.add(current.loci);
                this.mark(current, MarkerAction.Select);
            }
        }

        protected mark(current: Loci<ModelLoci>, action: MarkerAction) {
            for (let p of this.providers) p(current, action);
        }

        abstract apply(e: MarkEvent): void

        constructor(public ctx: PluginContext) {
            this.sel = ctx.helpers.structureSelection
        }
    }

    export class LociHighlightManager extends LociMarkManager<HighlightEvent> {
        private prev: Loci = { loci: EmptyLoci, repr: void 0 };

        apply(e: HighlightEvent) {
            const { current, modifiers } = e
            if (StructureElement.isLoci(current.loci)) {
                let loci: StructureElement.Loci = current.loci;
                if (modifiers && modifiers.shift) {
                    loci = this.sel.tryGetRange(loci) || loci;
                }

                this.mark(this.prev, MarkerAction.RemoveHighlight);
                const toHighlight = { loci, repr: current.repr };
                this.mark(toHighlight, MarkerAction.Highlight);
                this.prev = toHighlight;
            } else {
                if (!Loci.areEqual(this.prev, current)) {
                    this.mark(this.prev, MarkerAction.RemoveHighlight);
                    this.mark(current, MarkerAction.Highlight);
                    this.prev = current;
                }
            }
        }

        constructor(public ctx: PluginContext) {
            super(ctx)
            ctx.behaviors.interaction.highlight.subscribe(e => this.apply(e));
        }
    }

    export class LociSelectionManager extends LociMarkManager<ClickEvent> {
        toggleSel(current: Loci<ModelLoci>) {
            if (this.sel.has(current.loci)) {
                this.sel.remove(current.loci);
                this.mark(current, MarkerAction.Deselect);
            } else {
                this.sel.add(current.loci);
                this.mark(current, MarkerAction.Select);
            }
        }

        apply(e: ClickEvent) {
            const { current, buttons, modifiers } = e
            if (current.loci.kind === 'empty-loci') {
                if (modifiers.control && buttons === ButtonsType.Flag.Secondary) {
                    // clear the selection on Ctrl + Right-Click on empty
                    const sels = this.sel.clear();
                    for (const s of sels) this.mark({ loci: s }, MarkerAction.Deselect);
                }
            } else if (StructureElement.isLoci(current.loci)) {
                if (modifiers.control && buttons === ButtonsType.Flag.Secondary) {
                    // select only the current element on Ctrl + Right-Click
                    const old = this.sel.get(current.loci.structure);
                    this.mark({ loci: old }, MarkerAction.Deselect);
                    this.sel.set(current.loci);
                    this.mark(current, MarkerAction.Select);
                } else if (modifiers.control && buttons === ButtonsType.Flag.Primary) {
                    // toggle current element on Ctrl + Left-Click
                    this.toggleSel(current as Representation.Loci<StructureElement.Loci>);
                } else if (modifiers.shift && buttons === ButtonsType.Flag.Primary) {
                    // try to extend sequence on Shift + Left-Click
                    let loci: StructureElement.Loci = current.loci;
                    if (modifiers && modifiers.shift) {
                        loci = this.sel.tryGetRange(loci) || loci;
                    }
                    this.toggleSel({ loci, repr: current.repr });
                }
            } else {
                if (!ButtonsType.has(buttons, ButtonsType.Flag.Secondary)) return;
                for (let p of this.providers) p(current, MarkerAction.Toggle);
            }
        }

        constructor(public ctx: PluginContext) {
            super(ctx)
            ctx.behaviors.interaction.click.subscribe(e => this.apply(e));
        }
    }
}