/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PickingId } from '../../mol-geo/geometry/picking';
import { EmptyLoci } from '../../mol-model/loci';
import { Representation } from '../../mol-repr/representation';
import InputObserver, { ModifiersKeys, ButtonsType } from '../../mol-util/input/input-observer';
import { RxEventHelper } from '../../mol-util/rx-event-helper';

type Canvas3D = import('../canvas3d').Canvas3D
type HoverEvent = import('../canvas3d').Canvas3D.HoverEvent
type ClickEvent = import('../canvas3d').Canvas3D.ClickEvent

export class Canvas3dInteractionHelper {
    private ev = RxEventHelper.create();

    readonly events = {
        hover: this.ev<HoverEvent>(),
        click: this.ev<ClickEvent>(),
    };

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
    private button: ButtonsType.Flag = ButtonsType.create(0);
    private modifiers: ModifiersKeys = ModifiersKeys.None;

    private identify(isClick: boolean, t: number) {
        if (this.lastX !== this.cX || this.lastY !== this.cY) {
            this.id = this.canvasIdentify(this.cX, this.cY);
            this.lastX = this.cX;
            this.lastY = this.cY;
        }

        if (!this.id) return;

        if (isClick) {
            this.events.click.next({ current: this.getLoci(this.id), buttons: this.buttons, button: this.button, modifiers: this.modifiers });
            return;
        }

        if (!this.inside || this.currentIdentifyT !== t) {
            return;
        }

        const loci = this.getLoci(this.id);
        // only broadcast the latest hover
        if (!Representation.Loci.areEqual(this.prevLoci, loci)) {
            this.events.hover.next({ current: loci, buttons: this.buttons, button: this.button, modifiers: this.modifiers });
            this.prevLoci = loci;
        }
    }

    tick(t: number) {
        if (this.inside && t - this.prevT > 1000 / this.maxFps) {
            this.prevT = t;
            this.currentIdentifyT = t;
            this.identify(false, t);
        }
    }

    leave() {
        this.inside = false;
        if (this.prevLoci.loci !== EmptyLoci) {
            this.prevLoci = Representation.Loci.Empty;
            this.events.hover.next({ current: this.prevLoci, buttons: this.buttons, button: this.button, modifiers: this.modifiers });
        }
    }

    move(x: number, y: number, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys) {
        this.inside = true;
        this.buttons = buttons;
        this.button = button;
        this.modifiers = modifiers;
        this.cX = x;
        this.cY = y;
    }

    select(x: number, y: number, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys) {
        this.cX = x;
        this.cY = y;
        this.buttons = buttons;
        this.button = button;
        this.modifiers = modifiers;
        this.identify(true, 0);
    }

    modify(modifiers: ModifiersKeys) {
        if (this.prevLoci.loci === EmptyLoci || ModifiersKeys.areEqual(modifiers, this.modifiers)) return;
        this.modifiers = modifiers;
        this.events.hover.next({ current: this.prevLoci, buttons: this.buttons, button: this.button, modifiers: this.modifiers });
    }

    dispose() {
        this.ev.dispose();
    }

    constructor(private canvasIdentify: Canvas3D['identify'], private getLoci: Canvas3D['getLoci'], input: InputObserver, private maxFps: number = 30) {
        input.move.subscribe(({x, y, inside, buttons, button, modifiers }) => {
            if (!inside) return;
            this.move(x, y, buttons, button, modifiers);
        });

        input.leave.subscribe(() => {
            this.leave();
        });

        input.click.subscribe(({x, y, buttons, button, modifiers }) => {
            this.select(x, y, buttons, button, modifiers);
        });

        input.modifiers.subscribe(modifiers => this.modify(modifiers));
    }
}