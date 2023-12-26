/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ButtonsType, KeyCode, ModifiersKeys } from './input/input-observer';
import { camelCaseToWords, interpolate, stringToWords } from './string';

export { Binding };

interface Binding {
    triggers: Binding.Trigger[]
    action: string
    description: string
}

function Binding(triggers: Binding.Trigger[], action = '', description = '') {
    return Binding.create(triggers, action, description);
}

namespace Binding {
    export function create(triggers: Trigger[], action = '', description = ''): Binding {
        return { triggers, action, description };
    }

    export function isBinding(x: any): x is Binding {
        return !!x && Array.isArray(x.triggers) && typeof x.action === 'string';
    }

    export const Empty: Binding = { triggers: [], action: '', description: '' };
    export function isEmpty(binding: Binding) {
        return binding.triggers.length === 0 ||
            binding.triggers.every(t => t.buttons === undefined && t.modifiers === undefined && !t.code);
    }

    export function match(binding: Binding, buttons: ButtonsType, modifiers: ModifiersKeys) {
        return binding.triggers.some(t => Trigger.match(t, buttons, modifiers));
    }

    export function matchKey(binding: Binding, code: KeyCode, modifiers: ModifiersKeys, key: string) {
        return binding.triggers.some(t => Trigger.matchKey(t, code, modifiers, key));
    }

    export function formatTriggers(binding: Binding) {
        return binding.triggers.map(Trigger.format).join(' or ');
    }

    export function format(binding: Binding, name = '') {
        const help = binding.description || stringToWords(name);
        return interpolate(help, { triggers: '<i>' + formatTriggers(binding) + '</i>' });
    }

    export interface Trigger {
        buttons?: ButtonsType,
        modifiers?: ModifiersKeys
        code?: KeyCode
    }

    export function Trigger(buttons?: ButtonsType, modifiers?: ModifiersKeys) {
        return Trigger.create(buttons, modifiers);
    }

    export function TriggerKey(code?: KeyCode, modifiers?: ModifiersKeys) {
        return Trigger.create(undefined, modifiers, code);
    }

    export namespace Trigger {
        export function create(buttons?: ButtonsType, modifiers?: ModifiersKeys, code?: KeyCode): Trigger {
            return { buttons, modifiers, code };
        }
        export const Empty: Trigger = {};

        export function match(trigger: Trigger, buttons: ButtonsType, modifiers: ModifiersKeys): boolean {
            const { buttons: b, modifiers: m } = trigger;
            return b !== undefined &&
                (b === buttons || ButtonsType.has(b, buttons)) &&
                (!m || ModifiersKeys.areEqual(m, modifiers));
        }

        export function matchKey(trigger: Trigger, code: KeyCode, modifiers: ModifiersKeys, key: string): boolean {
            const { modifiers: m, code: c } = trigger;
            return c !== undefined &&
                (c === code || (
                    c.length === 1 &&
                    code.length === 4 &&
                    code.startsWith('Key') &&
                    !!key && key.length === 1 &&
                    key.toUpperCase() === c.toUpperCase()
                )) &&
                (!m || ModifiersKeys.areEqual(m, modifiers));
        }

        export function format(trigger: Trigger) {
            const s: string[] = [];
            const b = formatButtons(trigger.buttons, trigger.code);
            if (b) s.push(b);
            const c = formatCode(trigger.code);
            if (c) s.push(c);
            const m = formatModifiers(trigger.modifiers);
            if (m) s.push(m);
            return s.join(' + ');
        }
    }
}

const B = ButtonsType;

function formatButtons(buttons?: ButtonsType, code?: KeyCode) {
    const s: string[] = [];
    if (buttons === undefined && !code) {
        s.push('any mouse button');
    } else if (buttons === 0) {
        s.push('mouse hover');
    } else if (buttons !== undefined) {
        if (B.has(buttons, B.Flag.Primary)) s.push('left mouse button');
        if (B.has(buttons, B.Flag.Secondary)) s.push('right mouse button');
        if (B.has(buttons, B.Flag.Auxilary)) s.push('wheel/middle mouse button');
        if (B.has(buttons, B.Flag.Forth)) s.push('three fingers');
    }
    return s.join(' + ');
}

function formatModifiers(modifiers?: ModifiersKeys, verbose?: boolean) {
    const s: string[] = [];
    if (modifiers) {
        if (modifiers.alt) s.push('alt key');
        if (modifiers.control) s.push('control key');
        if (modifiers.meta) s.push('meta/command key');
        if (modifiers.shift) s.push('shift key');

        if (verbose && s.length === 0) s.push('no key');
    } else {
        if (verbose) s.push('any key');
    }
    return s.join(' + ');
}

function formatCode(code?: KeyCode) {
    if (code?.startsWith('Key')) code = code.substring(3);
    return code && camelCaseToWords(code).toLowerCase();
}