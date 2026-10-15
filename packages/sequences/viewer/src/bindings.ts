/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A behavior's bindings in Mol*'s terms (seqspec:viewer#behaviors): its params stay Seq*'s own (`Binding`,
// `KeyBinding`: plain data, matched allocation-free); for a help panel they are listed, with the provider's
// descriptions, as Mol*'s `Binding` (`Binding.format` renders one). Edge-rate.

import type { Binding as MolBinding } from '@molstar/core/util/binding';
import type { ButtonsType, ModifiersKeys } from '@molstar/core/util/input/input-observer';
import type { KeyTrigger, Trigger } from './input.js';

/** A provider's description of a binding param, by its path in the params (`"bindings.menu"`). */
export interface BindingDescription {
  readonly action: string;
  /** Mol*'s template: `${triggers}` is replaced by the triggers. */
  readonly description: string;
}
export type BindingDescriptions = { readonly [paramPath: string]: BindingDescription };
export interface SeqBindingInfo {
  readonly behavior: string;
  readonly param: string;
  readonly binding: MolBinding;
}

/** Mol*'s button flags by Seq*'s button name. */
const BUTTONS = { primary: 1, secondary: 2, middle: 4 } as const;

function trigger(t: Trigger | KeyTrigger): MolBinding.Trigger {
  const modifiers: ModifiersKeys = { shift: !!t.shift, control: !!t.ctrl, alt: !!t.alt, meta: !!t.meta };
  if ('key' in t || 'code' in t) return { modifiers, code: t.code ?? t.key };
  const b = (t as Trigger).button;
  return { modifiers, buttons: b === undefined ? undefined : (BUTTONS[b] as ButtonsType) };
}

/** The value at a dotted path of `params`; undefined: none. */
function at(params: object, path: string): unknown {
  let v: unknown = params;
  for (const k of path.split('.'))
    v = v !== null && typeof v === 'object' ? (v as Record<string, unknown>)[k] : undefined;
  return v;
}

/** One behavior's described bindings from its current params: a list of triggers, or one (`dominantAxis`); none: skipped. */
export function describeBindings(behavior: string, params: object, described: BindingDescriptions): SeqBindingInfo[] {
  const out: SeqBindingInfo[] = [];
  for (const [param, d] of Object.entries(described)) {
    const v = at(params, param);
    if (v === null || v === undefined || typeof v !== 'object') continue;
    const list = (Array.isArray(v) ? v : [v]) as readonly (Trigger | KeyTrigger)[];
    out.push({
      behavior,
      param,
      binding: { triggers: list.map(trigger), action: d.action, description: d.description },
    });
  }
  return out;
}
