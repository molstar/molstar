/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Sphere3D } from '../../../mol-math/geometry';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Loci } from '../../../mol-model/loci';
import { Structure, StructureElement } from '../../../mol-model/structure';
import { PluginContext } from '../../../mol-plugin/context';
import { StateSelection } from '../../../mol-state';
import { lociLabel } from '../../../mol-theme/label';
import { arrayRemoveAtInPlace } from '../../../mol-util/array';
import { StatefulPluginComponent } from '../../component';
import { PluginStateObject } from '../../objects';
import { getLociRange } from './selection';

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
    };

    readonly behaviors = {
        current: this.ev.behavior<FocusEntry | undefined>(void 0)
    };

    get current() { return this.state.current; }
    get history() { return this.state.history; }

    /** Last added or removed loci */
    private referenceLoci: StructureElement.Loci | undefined;

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
            this.events.historyUpdated.next(void 0);
            return;
        }

        this.state.history.unshift(entry);
        if (this.state.history.length > HISTORY_CAPACITY) this.state.history.pop();

        this.events.historyUpdated.next(void 0);
    }

    set(entry: FocusEntry) {
        this.tryAddHistory(entry);
        if (!this.state.current || !StructureElement.Loci.areEqual(this.state.current.loci, entry.loci)) {
            this.state.current = entry;
            this.behaviors.current.next(entry);
        }
    }

    tryGetRange(loci: Loci): StructureElement.Loci | undefined {
        return getLociRange(this.referenceLoci, loci);
    }

    setFromLoci(anyLoci: Loci) {
        const loci = Loci.normalize(anyLoci);
        if (!StructureElement.Loci.is(loci) || StructureElement.Loci.isEmpty(loci)) {
            this.clear();
            return;
        }

        this.set({ loci, label: lociLabel(loci, { reverse: true, hidePrefix: true, htmlStyling: false }) });
        this.referenceLoci = loci;
    }

    addFromLoci(anyLoci: Loci) {
        const union = this.state.current && StructureElement.Loci.is(anyLoci) && anyLoci.structure === this.state.current.loci.structure
            ? StructureElement.Loci.union(anyLoci, this.state.current.loci)
            : anyLoci;
        this.setFromLoci(union);
        const refLoci = Loci.normalize(anyLoci);
        this.referenceLoci = StructureElement.Loci.is(refLoci) ? refLoci : undefined;
    }

    toggleFromLoci(anyLoci: Loci) {
        const { kind, loci } = toggleLoci(this.state.current?.loci, anyLoci);
        this.setFromLoci(loci);
        const refLoci = Loci.normalize(anyLoci);
        this.referenceLoci = StructureElement.Loci.is(refLoci) && kind !== 'subtract' ? refLoci : undefined;
    }

    extendFromLoci(anyLoci: Loci) {
        const range = this.tryGetRange(anyLoci) ?? anyLoci;
        this.toggleFromLoci(range);

    }

    clear() {
        if (this.state.current) {
            this.state.current = undefined;
            this.behaviors.current.next(void 0);
        }
        this.referenceLoci = undefined;
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

    constructor(private plugin: PluginContext) {
        super({ history: [] });

        plugin.state.data.events.object.removed.subscribe(({ obj }) => {
            if (!PluginStateObject.Molecule.Structure.is(obj)) return;

            if (this.current?.loci.structure === obj.data) {
                this.clear();
            }

            const keep: FocusEntry[] = [];
            for (const e of this.history) {
                if (e.loci.structure === obj.data) keep.push(e);
            }
            if (keep.length !== this.history.length) {
                this.history.length = 0;
                this.history.push(...keep);
                this.events.historyUpdated.next(void 0);
            }
        });

        const sphere = Sphere3D();

        plugin.state.data.events.object.updated.subscribe(({ oldData, obj, action }) => {
            if (!PluginStateObject.Molecule.Structure.is(obj)) return;
            // structure NOT changed, keep everything as is; fixes #123
            if (oldData === obj.data) return;

            // structure changed (e.g. coordinates), try to remap and re-focus
            if (action === 'in-place') {
                const current = this.state.current;
                const structure = obj.data as Structure;

                if (current && current.loci.structure === oldData) {
                    const loci = StructureElement.Loci.remap(current.loci, structure);
                    this.state.current = { ...current, loci };
                    this.behaviors.current.next(this.state.current);

                    Loci.getBoundingSphere(loci, sphere);
                    const camera = this.plugin.canvas3d?.camera!;
                    const d = camera.getTargetDistance(sphere.radius + 4); // default extraRadius
                    if (Vec3.distance(camera.target, sphere.center) > sphere.radius ||
                        d > camera.viewport.height / camera.zoom
                    ) {
                        this.plugin.managers.camera.focusSphere(sphere, { durationMs: 0 });
                    }
                }

                // TODO remap history as well
            }
        });
    }
}

/** Return union of `currentLoci` and `newLoci`; or subtract `newLoci` from `currentLoci` if `newLoci` is a subset of `currentLoci`. */
function toggleLoci(currentLoci: StructureElement.Loci | undefined, newLoci: Loci) {
    if (currentLoci && StructureElement.Loci.is(newLoci) && newLoci.structure === currentLoci.structure) {
        if (StructureElement.Loci.isSubset(currentLoci, newLoci)) {
            return { kind: 'subtract', loci: StructureElement.Loci.subtract(currentLoci, newLoci) };
        } else {
            return { kind: 'add', loci: StructureElement.Loci.union(newLoci, currentLoci) };
        }
    } else {
        return { kind: 'new', loci: newLoci };
    }
}
