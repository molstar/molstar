/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PickingId } from '../../mol-geo/geometry/picking';
import { Representation } from '../../mol-repr/representation';
import { InputObserver, ModifiersKeys, ButtonsType } from '../../mol-util/input/input-observer';
import { RxEventHelper } from '../../mol-util/rx-event-helper';
import { Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { Camera } from '../camera';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Bond } from '../../mol-model/structure';
import { TrackballControls } from '../controls/trackball';

type Canvas3D = import('../canvas3d').Canvas3D
type HoverEvent = import('../canvas3d').Canvas3D.HoverEvent
type DragEvent = import('../canvas3d').Canvas3D.DragEvent
type ClickEvent = import('../canvas3d').Canvas3D.ClickEvent

enum InputEvent { Move, Click, Drag }

const tmpPosA = Vec3();
const tmpPos = Vec3();
const tmpNorm = Vec3();

export const Canvas3dInteractionHelperParams = {
    maxFps: PD.Numeric(30, { min: 10, max: 60, step: 10 }),
    preferAtomPixelPadding: PD.Numeric(3, { min: 0, max: 20, step: 1 }, { description: 'Number of extra pixels at which to prefer atoms over bonds.' }),
};
export type Canvas3dInteractionHelperParams = typeof Canvas3dInteractionHelperParams
export type Canvas3dInteractionHelperProps = PD.Values<Canvas3dInteractionHelperParams>

export class Canvas3dInteractionHelper {
    private ev = RxEventHelper.create();

    readonly events = {
        hover: this.ev<HoverEvent>(),
        drag: this.ev<DragEvent>(),
        click: this.ev<ClickEvent>(),
    };

    private startX = -1;
    private startY = -1;
    private endX = -1;
    private endY = -1;

    private id: PickingId | undefined = void 0;
    private position: Vec3 | undefined = void 0;

    private currentIdentifyT = 0;
    private isInteracting = false;

    private prevLoci: Representation.Loci = Representation.Loci.Empty;
    private prevT = 0;

    private inside = false;

    private buttons: ButtonsType = ButtonsType.create(0);
    private button: ButtonsType.Flag = ButtonsType.create(0);
    private modifiers: ModifiersKeys = ModifiersKeys.None;

    readonly props: Canvas3dInteractionHelperProps;

    setProps(props: Partial<Canvas3dInteractionHelperProps>) {
        Object.assign(this.props, props);
    }

    private identify(e: InputEvent, t: number) {
        const xyChanged = this.startX !== this.endX || this.startY !== this.endY || (this.input.pointerLock && !this.controls.isMoving);

        if (e === InputEvent.Drag) {
            if (xyChanged && !this.outsideViewport(this.startX, this.startY)) {
                this.events.drag.next({ current: this.prevLoci, buttons: this.buttons, button: this.button, modifiers: this.modifiers, pageStart: Vec2.create(this.startX, this.startY), pageEnd: Vec2.create(this.endX, this.endY) });

                this.startX = this.endX;
                this.startY = this.endY;
            }
            return;
        }

        if (xyChanged) {
            const pickData = this.canvasIdentify(this.endX, this.endY);
            this.id = pickData?.id;
            this.position = pickData?.position;
            this.startX = this.endX;
            this.startY = this.endY;
        }

        if (e === InputEvent.Click) {
            const loci = this.getLoci(this.id, this.position);
            this.events.click.next({ current: loci, buttons: this.buttons, button: this.button, modifiers: this.modifiers, page: Vec2.create(this.endX, this.endY), position: this.position });
            this.prevLoci = loci;
            return;
        }

        if (!this.inside || this.currentIdentifyT !== t || !xyChanged || this.outsideViewport(this.endX, this.endY)) return;

        const loci = this.getLoci(this.id, this.position);
        this.events.hover.next({ current: loci, buttons: this.buttons, button: this.button, modifiers: this.modifiers, page: Vec2.create(this.endX, this.endY), position: this.position });
        this.prevLoci = loci;
    }

    tick(t: number) {
        if (this.inside && t - this.prevT > 1000 / this.props.maxFps) {
            this.prevT = t;
            this.currentIdentifyT = t;
            this.identify(this.isInteracting ? InputEvent.Drag : InputEvent.Move, t);
        }
    }

    private leave() {
        this.inside = false;
        if (!Representation.Loci.isEmpty(this.prevLoci)) {
            this.prevLoci = Representation.Loci.Empty;
            this.events.hover.next({ current: this.prevLoci, buttons: this.buttons, button: this.button, modifiers: this.modifiers });
        }
    }

    private move(x: number, y: number, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys) {
        this.inside = true;
        this.buttons = buttons;
        this.button = button;
        this.modifiers = modifiers;
        this.endX = x;
        this.endY = y;
    }

    private click(x: number, y: number, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys) {
        this.endX = x;
        this.endY = y;
        this.buttons = buttons;
        this.button = button;
        this.modifiers = modifiers;
        this.identify(InputEvent.Click, 0);
    }

    private drag(x: number, y: number, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys) {
        this.endX = x;
        this.endY = y;
        this.buttons = buttons;
        this.button = button;
        this.modifiers = modifiers;
        this.identify(InputEvent.Drag, 0);
    }

    private modify(modifiers: ModifiersKeys) {
        if (ModifiersKeys.areEqual(modifiers, this.modifiers)) return;
        this.modifiers = modifiers;
        this.events.hover.next({ current: this.prevLoci, buttons: this.buttons, button: this.button, modifiers: this.modifiers, page: Vec2.create(this.endX, this.endY), position: this.position });
    }

    private outsideViewport(x: number, y: number) {
        const { input, camera: { viewport } } = this;
        x *= input.pixelRatio;
        y *= input.pixelRatio;
        return (
            x > viewport.x + viewport.width ||
            input.height - y > viewport.y + viewport.height ||
            x < viewport.x ||
            input.height - y < viewport.y
        );
    }

    private getLoci(pickingId: PickingId | undefined, position: Vec3 | undefined) {
        const { repr, loci } = this.lociGetter(pickingId);
        if (position && repr && Bond.isLoci(loci) && loci.bonds.length === 2) {
            const { aUnit, aIndex } = loci.bonds[0];
            aUnit.conformation.position(aUnit.elements[aIndex], tmpPosA);
            Vec3.sub(tmpNorm, this.camera.state.position, this.camera.state.target);
            Vec3.projectPointOnPlane(tmpPos, position, tmpNorm, tmpPosA);
            const pixelSize = this.camera.getPixelSize(tmpPos);
            let radius = repr.theme.size.size(loci.bonds[0]) * (repr.props.sizeFactor ?? 1);
            if (repr.props.lineSizeAttenuation === false) {
                // divide by two to get radius
                radius *= pixelSize / 2;
            }
            radius += this.props.preferAtomPixelPadding * pixelSize;
            if (Vec3.distance(tmpPos, tmpPosA) < radius) {
                return { repr, loci: Bond.toFirstStructureElementLoci(loci) };
            }
        }
        return { repr, loci };
    }

    dispose() {
        this.ev.dispose();
    }

    constructor(private canvasIdentify: Canvas3D['identify'], private lociGetter: Canvas3D['getLoci'], private input: InputObserver, private camera: Camera, private controls: TrackballControls, props: Partial<Canvas3dInteractionHelperProps> = {}) {
        this.props = { ...PD.getDefaultValues(Canvas3dInteractionHelperParams), ...props };

        input.drag.subscribe(({ x, y, buttons, button, modifiers }) => {
            this.isInteracting = true;
            // console.log('drag');
            this.drag(x, y, buttons, button, modifiers);
        });

        input.move.subscribe(({ x, y, inside, buttons, button, modifiers, onElement }) => {
            if (!inside || this.isInteracting) return;
            if (!onElement) {
                this.leave();
                return;
            }
            // console.log('move');
            this.move(x, y, buttons, button, modifiers);
        });

        input.leave.subscribe(() => {
            // console.log('leave');
            this.leave();
        });

        input.click.subscribe(({ x, y, buttons, button, modifiers }) => {
            if (this.outsideViewport(x, y)) return;
            // console.log('click');
            this.click(x, y, buttons, button, modifiers);
        });

        input.interactionEnd.subscribe(() => {
            // console.log('interactionEnd');
            this.isInteracting = false;
        });

        input.modifiers.subscribe(modifiers => {
            // console.log('modifiers');
            this.modify(modifiers);
        });
    }
}