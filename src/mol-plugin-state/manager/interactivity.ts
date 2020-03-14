/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Loci as ModelLoci, EmptyLoci, EveryLoci, isEmptyLoci } from '../../mol-model/loci';
import { ModifiersKeys, ButtonsType } from '../../mol-util/input/input-observer';
import { Representation } from '../../mol-repr/representation';
import { StructureElement } from '../../mol-model/structure';
import { MarkerAction } from '../../mol-util/marker-action';
import { PluginContext } from '../../mol-plugin/context';
import { Structure } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StructureSelectionManager } from './structure/selection';
import { PluginComponent } from '../component';
import { shallowEqual } from '../../mol-util/object';

export { InteractivityManager }

interface InteractivityManagerState {
    props: PD.ValuesFor<InteractivityManager.Params>
}

class InteractivityManager extends PluginComponent<InteractivityManagerState> {
    readonly lociSelects: InteractivityManager.LociSelectManager;
    readonly lociHighlights: InteractivityManager.LociHighlightManager;

    private _props = PD.getDefaultValues(InteractivityManager.Params)

    readonly events = {
        propsUpdated: this.ev()
    };

    get props(): Readonly<InteractivityManagerState['props']> { return { ...this.state.props }; }

    setProps(props: Partial<InteractivityManager.Props>) {
        const old = this.props;
        const _new = { ...this.state.props, ...props };
        if (shallowEqual(old, _new)) return;

        this.updateState({ props: _new });
        this.lociSelects.setProps(_new);
        this.lociHighlights.setProps(_new);
        this.events.propsUpdated.next();
    }

    constructor(readonly ctx: PluginContext, props: Partial<InteractivityManager.Props> = {}) {
        super({ props: { ...PD.getDefaultValues(InteractivityManager.Params), ...props } });

        this.lociSelects = new InteractivityManager.LociSelectManager(ctx, this._props);
        this.lociHighlights = new InteractivityManager.LociHighlightManager(ctx, this._props);
    }
}

namespace InteractivityManager {
    export interface Loci<T extends ModelLoci = ModelLoci> { loci: T, repr?: Representation.Any }

    export namespace Loci {
        export function areEqual(a: Loci, b: Loci) {
            return a.repr === b.repr && ModelLoci.areEqual(a.loci, b.loci);
        }
        export const Empty: Loci = { loci: EmptyLoci };
    }

    export const Params = {
        granularity: PD.Select('residue', ModelLoci.GranularityOptions, { label: 'Picking', description: 'Controls if selections are expanded upon picking to whole residues, chains, structures, instances, or left as atoms and coarse elements' }),
    }
    export type Params = typeof Params
    export type Props = PD.Values<Params>

    export interface HoverEvent { current: Loci, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys }
    export interface ClickEvent { current: Loci, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys }

    export type LociMarkProvider = (loci: Loci, action: MarkerAction) => void

    export abstract class LociMarkManager {
        protected providers: LociMarkProvider[] = [];
        protected sel: StructureSelectionManager

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

        protected normalizedLoci(interactivityLoci: Loci, applyGranularity = true) {
            const { loci, repr } = interactivityLoci
            const granularity =  applyGranularity ? this.props.granularity : undefined
            return { loci: ModelLoci.normalize(loci, granularity), repr }
        }

        protected mark(current: Loci<ModelLoci>, action: MarkerAction) {
            for (let p of this.providers) p(current, action);
        }

        constructor(public readonly ctx: PluginContext, props: Partial<Props> = {}) {
            this.sel = ctx.managers.structure.selection
            this.setProps(props)
        }
    }

    //

    export class LociHighlightManager extends LociMarkManager {
        private prev: Loci[] = [];

        private isHighlighted(loci: Loci) {
            for (const p of this.prev) {
                if (Loci.areEqual(p, loci)) return true
            }
            return false
        }

        private addHighlight(loci: Loci) {
            this.mark(loci, MarkerAction.Highlight);
            this.prev.push(loci)
        }

        clearHighlights() {
            for (const p of this.prev) {
                this.mark(p, MarkerAction.RemoveHighlight);
            }
            this.prev.length = 0
        }

        highlight(current: Loci, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity)
            if (!this.isHighlighted(normalized)) {
                this.addHighlight(normalized)
            }
        }

        highlightOnly(current: Loci, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity)
            if (!this.isHighlighted(normalized)) {
                this.clearHighlights()
                this.addHighlight(normalized)
            }
        }

        highlightOnlyExtend(current: Loci, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity)
            if (StructureElement.Loci.is(normalized.loci)) {
                const loci = {
                    loci: this.sel.tryGetRange(normalized.loci) || normalized.loci,
                    repr: normalized.repr
                }
                if (!this.isHighlighted(loci)) {
                    this.clearHighlights()
                    this.addHighlight(loci)
                }
            }
        }
    }

    //

    export class LociSelectManager extends LociMarkManager {
        toggle(current: Loci<ModelLoci>, applyGranularity = true) {
            if (ModelLoci.isEmpty(current.loci)) return;

            const normalized = this.normalizedLoci(current, applyGranularity)
            if (StructureElement.Loci.is(normalized.loci)) {
                this.toggleSel(normalized);
            } else {
                super.mark(normalized, MarkerAction.Toggle);
            }
        }

        toggleExtend(current: Loci<ModelLoci>, applyGranularity = true) {
            if (ModelLoci.isEmpty(current.loci)) return;

            const normalized = this.normalizedLoci(current, applyGranularity)
            if (StructureElement.Loci.is(normalized.loci)) {
                const loci = this.sel.tryGetRange(normalized.loci) || normalized.loci;
                this.toggleSel({ loci, repr: normalized.repr });
            }
        }

        select(current: Loci<ModelLoci>, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity)
            if (StructureElement.Loci.is(normalized.loci)) {
                this.sel.modify('add', normalized.loci);
            }
            this.mark(normalized, MarkerAction.Select);
        }

        selectJoin(current: Loci<ModelLoci>, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity)
            if (StructureElement.Loci.is(normalized.loci)) {
                this.sel.modify('intersect', normalized.loci);
            }
            this.mark(normalized, MarkerAction.Select);
        }

        selectOnly(current: Loci<ModelLoci>, applyGranularity = true) {
            this.deselectAll()
            const normalized = this.normalizedLoci(current, applyGranularity)
            if (StructureElement.Loci.is(normalized.loci)) {
                this.sel.modify('set', normalized.loci);
            }
            this.mark(normalized, MarkerAction.Select);
        }

        deselect(current: Loci<ModelLoci>, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity)
            if (StructureElement.Loci.is(normalized.loci)) {
                this.sel.modify('remove', normalized.loci);
            }
            this.mark(normalized, MarkerAction.Deselect);
        }

        deselectAll() {
            this.sel.clear();
            this.mark({ loci: EveryLoci }, MarkerAction.Deselect);
        }

        deselectAllOnEmpty(current: Loci<ModelLoci>) {
            if (isEmptyLoci(current.loci)) this.deselectAll()
        }

        protected mark(current: Loci<ModelLoci>, action: MarkerAction.Select | MarkerAction.Deselect) {
            const { loci } = current
            if (StructureElement.Loci.is(loci)) {
                // do a full deselect/select for the current structure so visuals that are
                // marked with granularity unequal to 'element' and join/intersect operations
                // are handled properly
                super.mark({ loci: Structure.Loci(loci.structure) }, MarkerAction.Deselect)
                super.mark({ loci: this.sel.getLoci(loci.structure) }, MarkerAction.Select)
            } else {
                super.mark(current, action)
            }
        }

        private toggleSel(current: Loci<ModelLoci>) {
            if (this.sel.has(current.loci)) {
                this.sel.modify('remove', current.loci);
                this.mark(current, MarkerAction.Deselect);
            } else {
                this.sel.modify('add', current.loci);
                this.mark(current, MarkerAction.Select);
            }
        }
    }
}