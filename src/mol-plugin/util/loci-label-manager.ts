/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from 'mol-plugin/context';
import { Loci } from 'mol-model/loci';
import { Representation } from 'mol-repr/representation';

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

    private empty: any[] = [];
    private getInfo(loci: Loci, repr?: Representation<any>) {
        if (!loci || loci.kind === 'empty-loci') return this.empty;
        const info: LociLabelEntry[] = [];
        for (let p of this.providers) {
            const e = p(loci, repr);
            if (e) info.push(e);
        }
        return info;
    }

    constructor(public ctx: PluginContext) {
        ctx.behaviors.canvas.highlightLoci.subscribe(ev => ctx.events.labels.highlight.next({ entries: this.getInfo(ev.loci, ev.repr) }));
    }
}