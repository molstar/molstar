/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { EveryLoci, isEmptyLoci, Loci } from '../../mol-model/loci';
import { Structure, StructureElement } from '../../mol-model/structure';
import { PluginContext } from '../../mol-plugin/context';
import { Representation } from '../../mol-repr/representation';
import { ButtonsType, ModifiersKeys } from '../../mol-util/input/input-observer';
import { MarkerAction } from '../../mol-util/marker-action';
import { shallowEqual } from '../../mol-util/object';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StatefulPluginComponent } from '../component';
import { StructureSelectionManager } from './structure/selection';

export { InteractivityManager };

interface InteractivityManagerState {
    props: PD.ValuesFor<InteractivityManager.Params>
}

// TODO: make this customizable somewhere?
const DefaultInteractivityFocusOptions = {
    minRadius: 6,
    extraRadius: 6,
    durationMs: 250,
};

export type InteractivityFocusLociOptions = typeof DefaultInteractivityFocusOptions

class InteractivityManager extends StatefulPluginComponent<InteractivityManagerState> {
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

    constructor(readonly plugin: PluginContext, props: Partial<InteractivityManager.Props> = {}) {
        super({ props: { ...PD.getDefaultValues(InteractivityManager.Params), ...props } });

        this.lociSelects = new InteractivityManager.LociSelectManager(plugin, this._props);
        this.lociHighlights = new InteractivityManager.LociHighlightManager(plugin, this._props);
    }
}

namespace InteractivityManager {
    export const Params = {
        granularity: PD.Select('residue', Loci.GranularityOptions, { label: 'Picking Level', description: 'Controls if selections are expanded upon picking to whole residues, chains, structures, instances, or left as atoms and coarse elements' }),
    };
    export type Params = typeof Params
    export type Props = PD.Values<Params>

    export interface HoverEvent { current: Representation.Loci, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys }
    export interface ClickEvent { current: Representation.Loci, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys }

    export type LociMarkProvider = (loci: Representation.Loci, action: MarkerAction) => void

    export abstract class LociMarkManager {
        protected providers: LociMarkProvider[] = [];
        protected sel: StructureSelectionManager

        readonly props: Readonly<Props> = PD.getDefaultValues(Params)

        setProps(props: Partial<Props>) {
            Object.assign(this.props, props);
        }

        addProvider(provider: LociMarkProvider) {
            this.providers.push(provider);
        }

        removeProvider(provider: LociMarkProvider) {
            this.providers = this.providers.filter(p => p !== provider);
            // TODO clear, then re-apply remaining providers
        }

        protected normalizedLoci(reprLoci: Representation.Loci, applyGranularity = true) {
            const { loci, repr } = reprLoci;
            const granularity =  applyGranularity ? this.props.granularity : undefined;
            return { loci: Loci.normalize(loci, granularity), repr };
        }

        protected mark(current: Representation.Loci, action: MarkerAction) {
            for (let p of this.providers) p(current, action);
        }

        constructor(public readonly ctx: PluginContext, props: Partial<Props> = {}) {
            this.sel = ctx.managers.structure.selection;
            this.setProps(props);
        }
    }

    //

    export class LociHighlightManager extends LociMarkManager {
        private prev: Representation.Loci[] = [];

        private isHighlighted(loci: Representation.Loci) {
            for (const p of this.prev) {
                if (Representation.Loci.areEqual(p, loci)) return true;
            }
            return false;
        }

        private addHighlight(loci: Representation.Loci) {
            this.mark(loci, MarkerAction.Highlight);
            this.prev.push(loci);
        }

        clearHighlights = () => {
            for (const p of this.prev) {
                this.mark(p, MarkerAction.RemoveHighlight);
            }
            this.prev.length = 0;
        }

        highlight(current: Representation.Loci, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity);
            if (!this.isHighlighted(normalized)) {
                this.addHighlight(normalized);
            }
        }

        highlightOnly(current: Representation.Loci, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity);
            if (!this.isHighlighted(normalized)) {
                this.clearHighlights();
                this.addHighlight(normalized);
            }
        }

        highlightOnlyExtend(current: Representation.Loci, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity);
            if (StructureElement.Loci.is(normalized.loci)) {
                const loci = {
                    loci: this.sel.tryGetRange(normalized.loci) || normalized.loci,
                    repr: normalized.repr
                };
                if (!this.isHighlighted(loci)) {
                    this.clearHighlights();
                    this.addHighlight(loci);
                }
            }
        }
    }

    //

    export class LociSelectManager extends LociMarkManager {
        toggle(current: Representation.Loci, applyGranularity = true) {
            if (Loci.isEmpty(current.loci)) return;

            const normalized = this.normalizedLoci(current, applyGranularity);
            if (StructureElement.Loci.is(normalized.loci)) {
                this.toggleSel(normalized);
            } else {
                super.mark(normalized, MarkerAction.Toggle);
            }
        }

        toggleExtend(current: Representation.Loci, applyGranularity = true) {
            if (Loci.isEmpty(current.loci)) return;

            const normalized = this.normalizedLoci(current, applyGranularity);
            if (StructureElement.Loci.is(normalized.loci)) {
                const loci = this.sel.tryGetRange(normalized.loci) || normalized.loci;
                this.toggleSel({ loci, repr: normalized.repr });
            }
        }

        select(current: Representation.Loci, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity);
            if (StructureElement.Loci.is(normalized.loci)) {
                this.sel.modify('add', normalized.loci);
            }
            this.mark(normalized, MarkerAction.Select);
        }

        selectJoin(current: Representation.Loci, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity);
            if (StructureElement.Loci.is(normalized.loci)) {
                this.sel.modify('intersect', normalized.loci);
            }
            this.mark(normalized, MarkerAction.Select);
        }

        selectOnly(current: Representation.Loci, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity);
            if (StructureElement.Loci.is(normalized.loci)) {
                // only deselect for the structure of the given loci
                this.deselect({ loci: Structure.toStructureElementLoci(normalized.loci.structure), repr: normalized.repr }, false);
                this.sel.modify('set', normalized.loci);
            }
            this.mark(normalized, MarkerAction.Select);
        }

        deselect(current: Representation.Loci, applyGranularity = true) {
            const normalized = this.normalizedLoci(current, applyGranularity);
            if (StructureElement.Loci.is(normalized.loci)) {
                this.sel.modify('remove', normalized.loci);
            }
            this.mark(normalized, MarkerAction.Deselect);
        }

        deselectAll() {
            this.sel.clear();
            this.mark({ loci: EveryLoci }, MarkerAction.Deselect);
        }

        deselectAllOnEmpty(current: Representation.Loci) {
            if (isEmptyLoci(current.loci)) this.deselectAll();
        }

        protected mark(current: Representation.Loci, action: MarkerAction.Select | MarkerAction.Deselect) {
            const { loci } = current;
            if (StructureElement.Loci.is(loci)) {
                // do a full deselect/select for the current structure so visuals that are
                // marked with granularity unequal to 'element' and join/intersect operations
                // are handled properly
                super.mark({ loci: Structure.Loci(loci.structure) }, MarkerAction.Deselect);
                super.mark({ loci: this.sel.getLoci(loci.structure) }, MarkerAction.Select);
            } else {
                super.mark(current, action);
            }
        }

        private toggleSel(current: Representation.Loci) {
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