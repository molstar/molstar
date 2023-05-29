/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../mol-plugin/context';
import { Loci } from '../../mol-model/loci';
import { Representation } from '../../mol-repr/representation';
import { MarkerAction } from '../../mol-util/marker-action';
import { arrayRemoveAtInPlace } from '../../mol-util/array';

// any represents React element. For compatibility to including the type
export type LociLabel = string | any
export type LociLabelProvider = {
    label: (loci: Loci, repr?: Representation<any>) => LociLabel | undefined
    group?: (entry: LociLabel) => string
    /** Labels from providers with higher priority are shown first */
    priority?: number
}

export class LociLabelManager {
    providers: LociLabelProvider[] = [];

    clearProviders() {
        this.providers = [];
        this.isDirty = true;
        this.showLabels();
    }

    addProvider(provider: LociLabelProvider) {
        this.providers.push(provider);
        this.providers.sort((a, b) => (b.priority || 0) - (a.priority || 0));
        this.isDirty = true;
        this.showLabels();
    }

    removeProvider(provider: LociLabelProvider) {
        this.providers = this.providers.filter(p => p !== provider);
        this.isDirty = true;
        this.showLabels();
    }

    private locis: Representation.Loci[] = [];

    private mark(loci: Representation.Loci, action: MarkerAction) {
        const idx = this.locis.findIndex(l => Representation.Loci.areEqual(loci, l));
        if (idx === -1 && action === MarkerAction.Highlight) {
            this.locis.push(loci);
            this.isDirty = true;
        } else if (idx !== -1 && action === MarkerAction.RemoveHighlight) {
            arrayRemoveAtInPlace(this.locis, idx);
            this.isDirty = true;
        }
    }

    private isDirty = false;
    private labels: LociLabel[] = [];
    private groupedLabels = new Map<string, LociLabel[]>();

    private showLabels() {
        this.ctx.behaviors.labels.highlight.next({ labels: this.getLabels() });
    }

    private getLabels() {
        if (this.isDirty) {
            this.groupedLabels.clear();
            this.labels.length = 0;
            for (const provider of this.providers) {
                for (const loci of this.locis) {
                    if (Loci.isEmpty(loci.loci)) continue;
                    const label = provider.label(loci.loci, loci.repr);
                    if (label) {
                        const hash = provider.group ? provider.group(label) : label.toString();
                        const group = this.groupedLabels.get(hash);
                        if (group) group.push(label);
                        else this.groupedLabels.set(hash, [label]);
                    }
                }
            }
            this.labels.length = 0;
            this.groupedLabels.forEach((group, hash) => {
                const count = group.length;
                const entry = count > 1 && group[0] !== group[1]
                    ? hash : group[0];

                this.labels.push(count === 1 ? entry : `${entry} <small>|| \u00D7 ${count}</small>`);
            });
            this.isDirty = false;
        }
        return this.labels;
    }

    constructor(public ctx: PluginContext) {
        ctx.managers.interactivity.lociHighlights.addProvider((loci, action, noRender) => {
            if (this.providers.length === 0) return;

            this.mark(loci, action);
            if (!noRender) this.showLabels();
        });
    }
}