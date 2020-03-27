/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StatefulPluginComponent } from '../../component';
import { PluginContext } from '../../../mol-plugin/context';
import { arrayRemoveAtInPlace } from '../../../mol-util/array';
import { StructureElement, Bond, Structure } from '../../../mol-model/structure';
import { Loci } from '../../../mol-model/loci';
import { lociLabel } from '../../../mol-theme/label';

export type FocusEntry = { label: string, loci: StructureElement.Loci, category?: string }

interface StructureFocusManagerState {
    current?: FocusEntry,
    history: FocusEntry[],
}

const HISTORY_CAPACITY = 12;

export class StructureFocusManager extends StatefulPluginComponent<StructureFocusManagerState> {

    readonly events = {
        changed: this.ev<undefined>(),
        historyUpdated: this.ev<undefined>()
    }

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
        this.tryAddHistory(entry)
        if (!this.state.current || !StructureElement.Loci.areEqual(this.state.current.loci, entry.loci)) {
            this.state.current = entry
            this.events.changed.next()
        }
    }

    setFromLoci(anyLoci: Loci) {
        let loci: StructureElement.Loci;
        if (StructureElement.Loci.is(anyLoci)) {
            loci = anyLoci;
        } else if (Bond.isLoci(anyLoci)) {
            loci = Bond.toStructureElementLoci(anyLoci);
        } else if (Structure.isLoci(anyLoci)) {
            loci = Structure.toStructureElementLoci(anyLoci.structure);
        } else {
            this.clear()
            return
        }

        if (StructureElement.Loci.isEmpty(loci)) {
            this.clear()
            return
        }

        this.set({ loci, label: lociLabel(loci, { reverse: true, hidePrefix: true, htmlStyling: false }) })
    }

    clear() {
        if (this.state.current) {
            this.state.current = undefined
            this.events.changed.next()
        }
    }

    // this.subscribeObservable(this.plugin.events.state.object.removed, o => {
    //     if (!PluginStateObject.Molecule.Structure.is(o.obj) || !StructureElement.Loci.is(lastLoci)) return;
    //     if (lastLoci.structure === o.obj.data) {
    //         lastLoci = EmptyLoci;
    //     }
    // });

    // this.subscribeObservable(this.plugin.events.state.object.updated, o => {
    //     if (!PluginStateObject.Molecule.Structure.is(o.oldObj) || !StructureElement.Loci.is(lastLoci)) return;
    //     if (lastLoci.structure === o.oldObj.data) {
    //         lastLoci = EmptyLoci;
    //     }
    // });

    constructor(plugin: PluginContext) {
        super({ history: [] });

        // plugin.state.data.events.object.removed.subscribe(e => this.onRemove(e.ref));
        // plugin.state.data.events.object.updated.subscribe(e => this.onUpdate(e.ref, e.oldObj, e.obj));
    }
}