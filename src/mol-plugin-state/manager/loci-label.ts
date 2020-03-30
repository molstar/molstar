/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../mol-plugin/context';
import { Loci } from '../../mol-model/loci';
import { Representation } from '../../mol-repr/representation';
import { MarkerAction } from '../../mol-util/marker-action';
import { arrayRemoveAtInPlace } from '../../mol-util/array';

export type LociLabelEntry = JSX.Element | string
export type LociLabelProvider = (info: Loci, repr?: Representation<any>) => LociLabelEntry | undefined

export class LociLabelManager {
    providers: LociLabelProvider[] = [];

    addProvider(provider: LociLabelProvider) {
        this.providers.push(provider);
        this.isDirty = true
        this.showLabels()
    }

    removeProvider(provider: LociLabelProvider) {
        this.providers = this.providers.filter(p => p !== provider);
        this.isDirty = true
        this.showLabels()
    }

    private locis: Representation.Loci[] = []

    private mark(loci: Representation.Loci, action: MarkerAction) {
        const idx = this.locis.findIndex(l => Representation.Loci.areEqual(loci, l))
        if (idx === -1 && action === MarkerAction.Highlight) {
            this.locis.push(loci)
            this.isDirty = true
        } else if(idx !== -1 && action === MarkerAction.RemoveHighlight) {
            arrayRemoveAtInPlace(this.locis, idx)
            this.isDirty = true
        }
    }

    private isDirty = false
    private entries: LociLabelEntry[] = []
    private entriesCounts = new Map<LociLabelEntry, number>()

    private showLabels() {
        this.ctx.behaviors.labels.highlight.next({ entries: this.getEntries() })
    }

    private getEntries() {
        if (this.isDirty) {
            this.entriesCounts.clear()
            this.entries.length = 0
            for (const provider of this.providers) {
                for (const loci of this.locis) {
                    if (Loci.isEmpty(loci.loci)) continue
                    const entry = provider(loci.loci, loci.repr)
                    if (entry) {
                        const count = this.entriesCounts.get(entry) || 0
                        this.entriesCounts.set(entry, count + 1)
                    }
                }
            }
            this.entries.length = 0
            this.entriesCounts.forEach((count, entry) => {
                this.entries.push(count === 1 ? entry : `${entry} (\u00D7 ${count})`)
            })
            this.isDirty = false
        }
        return this.entries
    }

    constructor(public ctx: PluginContext) {
        ctx.managers.interactivity.lociHighlights.addProvider((loci, action) => {
            this.mark(loci, action)
            this.showLabels()
        })
    }
}