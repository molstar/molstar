/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StatefulPluginComponent } from '../../component';
import { PluginContext } from '../../../mol-plugin/context';
import { arrayRemoveAtInPlace } from '../../../mol-util/array';
import { StructureElement, Structure } from '../../../mol-model/structure';
import { Loci } from '../../../mol-model/loci';
import { lociLabel } from '../../../mol-theme/label';
import { PluginStateObject } from '../../objects';
import { StateSelection } from '../../../mol-state';

export type FocusEntry = {
    label: string
    loci: StructureElement.Loci
    category?: string
}

export interface StructureFocusSnapshot {
    current?: {
        label: string
        ref: string
        bundle: StructureElement.Bundle
        category?: string
    }
}

interface StructureFocusManagerState {
    current?: FocusEntry
    history: FocusEntry[]
}

const HISTORY_CAPACITY = 8;

export class StructureFocusManager extends StatefulPluginComponent<StructureFocusManagerState> {
    readonly events = {
        historyUpdated: this.ev<undefined>()
    }

    readonly behaviors = {
        current: this.ev.behavior<FocusEntry | undefined>(void 0)
    };

    get current() { return this.state.current; }
    get history() { return this.state.history; }

    private tryAddHistory(entry: FocusEntry) {
        if (StructureElement.Loci.isEmpty(entry.loci)) return;

        let idx = 0, existingEntry: FocusEntry | undefined = void 0;
        for (const e of this.state.history) {
            if (StructureElement.Loci.areEqual(e.loci, entry.loci)) {
                existingEntry = e;
                break;
            }
            idx++;
        }

        if (existingEntry) {
            // move to top, use new
            arrayRemoveAtInPlace(this.state.history, idx);
            this.state.history.unshift(entry);
            this.events.historyUpdated.next();
            return;
        }

        this.state.history.unshift(entry);
        if (this.state.history.length > HISTORY_CAPACITY) this.state.history.pop();

        this.events.historyUpdated.next();
    }

    set(entry: FocusEntry) {
        this.tryAddHistory(entry);
        if (!this.state.current || !StructureElement.Loci.areEqual(this.state.current.loci, entry.loci)) {
            this.state.current = entry;
            this.behaviors.current.next(entry);
        }
    }

    setFromLoci(anyLoci: Loci) {
        const loci = Loci.normalize(anyLoci);
        if (!StructureElement.Loci.is(loci) || StructureElement.Loci.isEmpty(loci)) {
            this.clear();
            return;
        }

        this.set({ loci, label: lociLabel(loci, { reverse: true, hidePrefix: true, htmlStyling: false }) });
    }

    addFromLoci(anyLoci: Loci) {
        const union = this.state.current && StructureElement.Loci.is(anyLoci) && anyLoci.structure === this.state.current.loci.structure
            ? StructureElement.Loci.union(anyLoci, this.state.current.loci)
            : anyLoci;
        this.setFromLoci(union);
    }

    clear() {
        if (this.state.current) {
            this.state.current = undefined;
            this.behaviors.current.next(void 0);
        }
    }

    getSnapshot(): StructureFocusSnapshot {
        if (!this.current) return {};

        const node = this.plugin.helpers.substructureParent.get(this.current.loci.structure);
        const ref = node?.transform.ref;
        if (!ref) return {};

        return {
            current: {
                label: this.current.label,
                ref,
                bundle: StructureElement.Bundle.fromLoci(this.current.loci),
                category: this.current.category
            }
        };
    }

    setSnapshot(snapshot: StructureFocusSnapshot) {
        if (!snapshot.current) {
            this.clear();
            return;
        }

        const { label, ref, bundle, category } = snapshot.current;
        const structure = this.plugin.state.data.select(StateSelection.Generators.byRef(ref))[0]?.obj?.data as Structure;
        if (!structure) return;

        const loci = StructureElement.Bundle.toLoci(bundle, structure);
        this.set({ label, loci, category });
    }

    // this.subscribeObservable(this.plugin.state.events.object.updated, o => {
    //     if (!PluginStateObject.Molecule.Structure.is(o.oldObj) || !StructureElement.Loci.is(lastLoci)) return;
    //     if (lastLoci.structure === o.oldObj.data) {
    //         lastLoci = EmptyLoci;
    //     }
    // });

    constructor(private plugin: PluginContext) {
        super({ history: [] });

        plugin.state.data.events.object.removed.subscribe(o => {
            if (!PluginStateObject.Molecule.Structure.is(o.obj)) return;

            if (this.current?.loci.structure === o.obj.data) {
                this.clear();
            }

            const keep: FocusEntry[] = [];
            for (const e of this.history) {
                if (e.loci.structure === o.obj.data) keep.push(e);
            }
            if (keep.length !== this.history.length) {
                this.history.length = 0;
                this.history.push(...keep);
                this.events.historyUpdated.next();
            }
        });
        // plugin.state.data.events.object.updated.subscribe(e => this.onUpdate(e.ref, e.oldObj, e.obj));
    }
}