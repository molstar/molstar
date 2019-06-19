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
import { StructureElement as SE, Structure } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginCommands } from '../command';
import { capitalize } from '../../mol-util/string';

export { Interactivity }

class Interactivity {
    readonly lociSelections: Interactivity.LociSelectionManager;
    readonly lociHighlights: Interactivity.LociHighlightManager;

    private _props = PD.getDefaultValues(Interactivity.Params)

    get props() { return { ...this._props } }
    setProps(props: Partial<Interactivity.Props>) {
        Object.assign(this._props, props)
        this.lociSelections.setProps(this._props)
        this.lociHighlights.setProps(this._props)
    }

    constructor(readonly ctx: PluginContext, props: Partial<Interactivity.Props> = {}) {
        Object.assign(this._props, props)

        this.lociSelections = new Interactivity.LociSelectionManager(ctx, this._props);
        this.lociHighlights = new Interactivity.LociHighlightManager(ctx, this._props);

        PluginCommands.Interactivity.SetProps.subscribe(ctx, e => this.setProps(e.props));
    }
}

namespace Interactivity {
    export interface Loci<T extends ModelLoci = ModelLoci> { loci: T, repr?: Representation.Any }

    export namespace Loci {
        export function areEqual(a: Loci, b: Loci) {
            return a.repr === b.repr && ModelLoci.areEqual(a.loci, b.loci);
        }
        export const Empty: Loci = { loci: EmptyLoci };
    }

    const LociExpansion = {
        'none': (loci: ModelLoci) => loci,
        'residue': (loci: ModelLoci) => SE.isLoci(loci) ? SE.Loci.extendToWholeResidues(loci) : loci,
        'chain': (loci: ModelLoci) => SE.isLoci(loci) ? SE.Loci.extendToWholeChains(loci) : loci,
        'structure': (loci: ModelLoci) => SE.isLoci(loci) ? Structure.Loci(loci.structure) : loci
    }
    type LociExpansion = keyof typeof LociExpansion
    const LociExpansionOptions = Object.keys(LociExpansion).map(n => [n, capitalize(n)]) as [LociExpansion, string][]

    export const Params = {
        lociExpansion: PD.Select('residue', LociExpansionOptions),
    }
    export type Props = PD.Values<typeof Params>

    export interface HighlightEvent { current: Loci, modifiers?: ModifiersKeys }
    export interface ClickEvent { current: Loci, buttons: ButtonsType, modifiers: ModifiersKeys }

    export type LociMarkProvider = (loci: Loci, action: MarkerAction) => void

    export abstract class LociMarkManager<MarkEvent extends any> {
        protected providers: LociMarkProvider[] = [];
        protected sel: StructureElementSelectionManager

        readonly props: Readonly<Props> = PD.getDefaultValues(Params)

        setProps(props: Partial<Props>) {
            Object.assign(this.props, props)
        }

        addProvider(provider: LociMarkProvider) {
            this.providers.push(provider);
        }

        removeProvider(provider: LociMarkProvider) {
            this.providers = this.providers.filter(p => p !== provider);
            // TODO clear, then re-apply remaining providers
        }

        expandLoci(loci: ModelLoci) {
            return LociExpansion[this.props.lociExpansion](loci)
        }

        protected mark(current: Loci<ModelLoci>, action: MarkerAction) {
            for (let p of this.providers) p(current, action);
        }

        abstract apply(e: MarkEvent): void

        constructor(public readonly ctx: PluginContext, props: Partial<Props> = {}) {
            this.sel = ctx.helpers.structureSelection
            this.setProps(props)
        }
    }

    export class LociHighlightManager extends LociMarkManager<HighlightEvent> {
        private prev: Loci = { loci: EmptyLoci, repr: void 0 };

        apply(e: HighlightEvent) {
            const { current, modifiers } = e
            const expanded: Loci<ModelLoci> = { loci: this.expandLoci(current.loci), repr: current.repr }
            if (StructureElement.isLoci(expanded.loci)) {
                let loci: StructureElement.Loci = expanded.loci;
                if (modifiers && modifiers.shift) {
                    loci = this.sel.tryGetRange(loci) || loci;
                }

                this.mark(this.prev, MarkerAction.RemoveHighlight);
                const toHighlight = { loci, repr: expanded.repr };
                this.mark(toHighlight, MarkerAction.Highlight);
                this.prev = toHighlight;
            } else {
                if (!Loci.areEqual(this.prev, expanded)) {
                    this.mark(this.prev, MarkerAction.RemoveHighlight);
                    this.mark(expanded, MarkerAction.Highlight);
                    this.prev = expanded;
                }
            }
        }

        constructor(ctx: PluginContext, props: Partial<Props> = {}) {
            super(ctx, props)
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
            const expanded: Loci<ModelLoci> = { loci: this.expandLoci(current.loci), repr: current.repr }
            if (expanded.loci.kind === 'empty-loci') {
                if (modifiers.control && buttons === ButtonsType.Flag.Secondary) {
                    // clear the selection on Ctrl + Right-Click on empty
                    const sels = this.sel.clear();
                    for (const s of sels) this.mark({ loci: s }, MarkerAction.Deselect);
                }
            } else if (StructureElement.isLoci(expanded.loci)) {
                if (modifiers.control && buttons === ButtonsType.Flag.Secondary) {
                    // select only the current element on Ctrl + Right-Click
                    const old = this.sel.get(expanded.loci.structure);
                    this.mark({ loci: old }, MarkerAction.Deselect);
                    this.sel.set(expanded.loci);
                    this.mark(expanded, MarkerAction.Select);
                } else if (modifiers.control && buttons === ButtonsType.Flag.Primary) {
                    // toggle current element on Ctrl + Left-Click
                    this.toggleSel(expanded as Representation.Loci<StructureElement.Loci>);
                } else if (modifiers.shift && buttons === ButtonsType.Flag.Primary) {
                    // try to extend sequence on Shift + Left-Click
                    let loci: StructureElement.Loci = expanded.loci;
                    if (modifiers && modifiers.shift) {
                        loci = this.sel.tryGetRange(loci) || loci;
                    }
                    this.toggleSel({ loci, repr: expanded.repr });
                }
            } else {
                if (!ButtonsType.has(buttons, ButtonsType.Flag.Secondary)) return;
                for (let p of this.providers) p(expanded, MarkerAction.Toggle);
            }
        }

        constructor(ctx: PluginContext, props: Partial<Props> = {}) {
            super(ctx, props)
            ctx.behaviors.interaction.click.subscribe(e => this.apply(e));
        }
    }
}