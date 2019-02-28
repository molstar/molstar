/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../context';
import { PickingId } from 'mol-geo/geometry/picking';
import { EmptyLoci } from 'mol-model/loci';
import { Representation } from 'mol-repr/representation';
import { ModifiersKeys, ButtonsType } from 'mol-util/input/input-observer';

export class Canvas3dIdentifyHelper {
    private cX = -1;
    private cY = -1;

    private lastX = -1;
    private lastY = -1;

    private id: PickingId | undefined = void 0;

    private currentIdentifyT = 0;

    private prevLoci: Representation.Loci = Representation.Loci.Empty;
    private prevT = 0;

    private inside = false;

    private buttons: ButtonsType = ButtonsType.create(0);
    private modifiers: ModifiersKeys = ModifiersKeys.None;

    private async identify(isClick: boolean, t: number) {
        if (this.lastX !== this.cX && this.lastY !== this.cY) {
            this.id = await this.ctx.canvas3d.identify(this.cX, this.cY);
            this.lastX = this.cX;
            this.lastY = this.cY;
        }

        if (!this.id) return;

        if (isClick) {
            this.ctx.events.canvas3d.click.next({ loci: this.ctx.canvas3d.getLoci(this.id), buttons: this.buttons, modifiers: this.modifiers });
            return;
        }

        // only highlight the latest
        if (!this.inside || this.currentIdentifyT !== t) {
            return;
        }

        const loci = this.ctx.canvas3d.getLoci(this.id);
        if (!Representation.Loci.areEqual(this.prevLoci, loci)) {
            this.ctx.events.canvas3d.highlight.next({ loci, modifiers: this.modifiers });
            this.prevLoci = loci;
        }
    }

    private animate: (t: number) => void = t => {
        if (!this.ctx.state.animation.isAnimating && this.inside && t - this.prevT > 1000 / this.maxFps) {
            this.prevT = t;
            this.currentIdentifyT = t;
            this.identify(false, t);
        }
        requestAnimationFrame(this.animate);
    }

    leave() {
        this.inside = false;
        if (this.prevLoci.loci !== EmptyLoci) {
            this.prevLoci = Representation.Loci.Empty;
            this.ctx.events.canvas3d.highlight.next({ loci: this.prevLoci });
            this.ctx.canvas3d.requestDraw(true);
        }
    }

    move(x: number, y: number, modifiers: ModifiersKeys) {
        this.inside = true;
        this.modifiers = modifiers;
        this.cX = x;
        this.cY = y;
    }

    select(x: number, y: number, buttons: ButtonsType, modifiers: ModifiersKeys) {
        this.cX = x;
        this.cY = y;
        this.buttons = buttons;
        this.modifiers = modifiers;
        this.identify(true, 0);
    }

    constructor(private ctx: PluginContext, private maxFps: number = 15) {
        this.animate(0);
    }
}