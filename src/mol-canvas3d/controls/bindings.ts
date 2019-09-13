/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ButtonsType, ModifiersKeys } from '../../mol-util/input/input-observer';

const B = ButtonsType
const M = ModifiersKeys

export interface Bindings {
    drag: {
        rotate: Bindings.Trigger
        rotateZ: Bindings.Trigger
        pan: Bindings.Trigger
        focus: Bindings.Trigger
        zoom: Bindings.Trigger
    },
    scroll: {
        focus: Bindings.Trigger
        zoom: Bindings.Trigger
    }
}

export namespace Bindings {
    export type Trigger = { buttons: ButtonsType, modifiers?: ModifiersKeys }

    export function match(trigger: Trigger, buttons: ButtonsType, modifiers: ModifiersKeys) {
        const { buttons: b, modifiers: m } = trigger
        return ButtonsType.has(b, buttons) &&
            (!m || ModifiersKeys.areEqual(m, modifiers))
    }

    export const Default: Bindings = {
        drag: {
            rotate: { buttons: B.Flag.Primary, modifiers: M.create() },
            rotateZ: { buttons: B.Flag.Primary, modifiers: M.create({ shift: true }) },
            pan: { buttons: B.Flag.Secondary, modifiers: M.create() },
            focus: { buttons: B.Flag.Forth, modifiers: M.create() },
            zoom: { buttons: B.Flag.Auxilary, modifiers: M.create() },
        },
        scroll: {
            focus: { buttons: B.Flag.Auxilary, modifiers: M.create({ shift: true }) },
            zoom: { buttons: B.Flag.Auxilary, modifiers: M.create() },
        }
    }
}