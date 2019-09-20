/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ButtonsType, ModifiersKeys } from './input/input-observer';
import { interpolate, stringToWords } from './string';

export { Binding }

interface Binding {
    trigger: Binding.Trigger
    description: string
}

function Binding(trigger: Binding.Trigger, description = '') {
    return Binding.create(trigger, description)
}

namespace Binding {
    export function create(trigger: Trigger, description = ''): Binding {
        return { trigger, description }
    }

    export const Empty: Binding = { trigger: {}, description: '' }
    export function isEmpty(binding: Binding) {
        return binding.trigger.buttons === undefined && binding.trigger.modifiers === undefined
    }

    export function match(binding: Binding, buttons: ButtonsType, modifiers: ModifiersKeys) {
        return Trigger.match(binding.trigger, buttons, modifiers)
    }

    export function format(binding: Binding, name = '') {
        const help = binding.description || stringToWords(name)
        return interpolate(help, { trigger: Trigger.format(binding.trigger) })
    }

    export interface Trigger {
        buttons?: ButtonsType,
        modifiers?: ModifiersKeys
    }

    export function Trigger(buttons?: ButtonsType, modifiers?: ModifiersKeys) {
        return Trigger.create(buttons, modifiers)
    }

    export namespace Trigger {
        export function create(buttons?: ButtonsType, modifiers?: ModifiersKeys): Trigger {
            return { buttons, modifiers }
        }
        export const Empty: Trigger = {}

        export function match(trigger: Trigger, buttons: ButtonsType, modifiers: ModifiersKeys): boolean {
            const { buttons: b, modifiers: m } = trigger
            return b !== undefined &&
                (b === buttons || ButtonsType.has(b, buttons)) &&
                (!m || ModifiersKeys.areEqual(m, modifiers))
        }

        export function format(trigger: Trigger) {
            const s: string[] = []
            const b = formatButtons(trigger.buttons)
            if (b) s.push(b)
            const m = formatModifiers(trigger.modifiers)
            if (m) s.push(m)
            return s.join(' + ')
        }
    }
}

const B = ButtonsType

function formatButtons(buttons?: ButtonsType) {
    const s: string[] = []
    if (buttons === undefined) {
        s.push('any button')
    } else if (buttons === 0) {
        s.push('no button')
    } else {
        if (B.has(buttons, B.Flag.Primary)) s.push('left button')
        if (B.has(buttons, B.Flag.Secondary)) s.push('right button')
        if (B.has(buttons, B.Flag.Auxilary)) s.push('wheel/middle button')
        if (B.has(buttons, B.Flag.Forth)) s.push('three fingers')
    }
    return s.join(' + ')
}

function formatModifiers(modifiers?: ModifiersKeys) {
    const s: string[] = []
    if (modifiers) {
        if (modifiers.alt) s.push('alt')
        if (modifiers.control) s.push('control')
        if (modifiers.meta) s.push('meta/command')
        if (modifiers.shift) s.push('shift')

        if (s.length === 0) s.push('no key')
    } else {
        s.push('any key')
    }
    return s.join(' + ')
}