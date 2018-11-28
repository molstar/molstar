/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../context';
import { PickingId } from 'mol-geo/geometry/picking';
import { EmptyLoci, Loci } from 'mol-model/loci';
import { Representation } from 'mol-repr/representation';

export class Canvas3dIdentifyHelper {
    private cX = -1;
    private cY = -1;

    private lastX = -1;
    private lastY = -1;

    private id: PickingId | undefined = void 0;

    private currentIdentifyT = 0;

    private prevLoci: { loci: Loci, repr?: Representation.Any } = { loci: EmptyLoci };
    private prevT = 0;

    private inside = false;

    private async identify(select: boolean, t: number) {
        if (this.lastX !== this.cX && this.lastY !== this.cY) {
            this.id = await this.ctx.canvas3d.identify(this.cX, this.cY);
            this.lastX = this.cX;
            this.lastY = this.cY;
        }

        if (!this.id) return;

        if (select) {
            this.ctx.behaviors.canvas.selectLoci.next(this.ctx.canvas3d.getLoci(this.id));
            return;
        }

        // only highlight the latest
        if (!this.inside || this.currentIdentifyT !== t) {
            return;
        }

        const loci = this.ctx.canvas3d.getLoci(this.id);
        if (loci.repr !== this.prevLoci.repr || !Loci.areEqual(loci.loci, this.prevLoci.loci)) {
            this.ctx.behaviors.canvas.highlightLoci.next(loci);
            this.prevLoci = loci;
        }
    }

    private animate: (t: number) => void = t => {
        if (this.inside && t - this.prevT > 1000 / this.maxFps) {
            this.prevT = t;
            this.currentIdentifyT = t;
            this.identify(false, t);
        }
        requestAnimationFrame(this.animate);
    }

    leave() {
        this.inside = false;
        if (this.prevLoci.loci !== EmptyLoci) {
            this.prevLoci = { loci: EmptyLoci };
            this.ctx.behaviors.canvas.highlightLoci.next(this.prevLoci);
            this.ctx.canvas3d.requestDraw(true);
        }
    }

    move(x: number, y: number) {
        this.inside = true;
        this.cX = x;
        this.cY = y;
    }

    select(x: number, y: number) {
        this.cX = x;
        this.cY = y;
        this.identify(true, 0);
    }

    constructor(private ctx: PluginContext, private maxFps: number = 15) {
        this.animate(0);
    }
}