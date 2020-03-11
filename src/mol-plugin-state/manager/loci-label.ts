/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../mol-plugin/context';
import { Loci } from '../../mol-model/loci';
import { Representation } from '../../mol-repr/representation';
import { MarkerAction } from '../../mol-util/marker-action';

export type LociLabelEntry = JSX.Element | string
export type LociLabelProvider = (info: Loci, repr?: Representation<any>) => LociLabelEntry | undefined

export class LociLabelManager {
    providers: LociLabelProvider[] = [];

    addProvider(provider: LociLabelProvider) {
        this.providers.push(provider);
    }

    removeProvider(provider: LociLabelProvider) {
        this.providers = this.providers.filter(p => p !== provider);
        // Event.Interactivity.Highlight.dispatch(this.ctx, []);
    }

    private empty: LociLabelEntry[] = [];
    private getInfo({ loci, repr }: Representation.Loci, action: MarkerAction) {
        if (Loci.isEmpty(loci) || action !== MarkerAction.Highlight) return this.empty;
        const info: LociLabelEntry[] = [];
        for (let p of this.providers) {
            const e = p(loci, repr);
            if (e) info.push(e);
        }
        return info;
    }

    constructor(public ctx: PluginContext) {
        ctx.managers.interactivity.lociHighlights.addProvider((loci, action) => ctx.behaviors.labels.highlight.next({ entries: this.getInfo(loci, action) }))
    }
}