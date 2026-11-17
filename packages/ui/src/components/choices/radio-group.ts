/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { KeyedSelector, type ReadonlyState } from '@molstar/core/reactive';
import { setAttr, watch } from '../../dom/index.js';

export interface RadioGroupOptions<K> {
  /** Each radio's key, in order: the value `active` holds when it is checked. */
  keys: readonly K[];
  /** The radios (`role="radio"`), in the order of `keys`. */
  buttons: readonly HTMLElement[];
  /** Whether the radio at an index is not available: never the fallback stop, skipped by the arrows. */
  disabled(index: number): boolean;
  /** The key in effect (none checked when it is no radio's key). */
  active: ReadonlyState<K | undefined>;
  /** Picks the radio at an index (an enabled one): the arrows, Home and End call it with the one they move to. */
  pick(index: number): void;
}

/**
 * A WAI-ARIA radio group over a row of buttons (APG radio group): marks the checked one (`aria-checked`) and keeps one
 * tab stop (a roving `tabindex`: the checked radio, else the first enabled), through one `KeyedSelector` over `active`,
 * so a change rewrites the two radios it concerns. Returns the group's `keydown` handler: the arrows move to the next or
 * previous enabled radio, wrapping, Home and End to the first and the last; the radio moved to takes the focus and is
 * picked. Call it in a view scope (it watches).
 */
export function radioGroup<K>(options: RadioGroupOptions<K>): (e: KeyboardEvent) => void {
  const { keys, buttons, active, disabled } = options;
  const n = keys.length;
  const enabled = (i: number) => !disabled(i);
  const selected = new KeyedSelector(active);
  const enabledKeys = new Set(keys.filter((_, i) => enabled(i)));
  const first = keys.findIndex((_, i) => enabled(i));
  const write = (el: HTMLElement, checked: boolean, stop: boolean) => {
    setAttr(el, 'aria-checked', String(checked));
    const tab = stop ? 0 : -1;
    if (el.tabIndex !== tab) el.tabIndex = tab;
  };
  keys.forEach((k, i) => {
    const el = buttons[i];
    const on = selected.key(k);
    if (i === first) {
      // The fallback stop: the first enabled radio holds it while no enabled radio is checked. Only it follows the
      // whole value; the others follow their own key.
      watch([on, active], () => write(el, on.value, on.value || !enabledKeys.has(active.value as K)));
    } else {
      const off = disabled(i);
      watch([on], () => write(el, on.value, on.value && !off));
    }
  });

  return (e: KeyboardEvent) => {
    if (e.defaultPrevented || e.altKey || e.ctrlKey || e.metaKey) return;
    const target = e.target as Element;
    const at = buttons.findIndex((b) => b === target || b.contains(target));
    if (at < 0) return;
    const step = (dir: number) => {
      for (let k = 1; k <= n; k++) {
        const j = (((at + dir * k) % n) + n) % n;
        if (enabled(j)) return j;
      }
      return -1;
    };
    let next: number;
    switch (e.key) {
      case 'ArrowRight':
      case 'ArrowDown':
        next = step(1);
        break;
      case 'ArrowLeft':
      case 'ArrowUp':
        next = step(-1);
        break;
      case 'Home':
        next = first;
        break;
      case 'End':
        next = n - 1;
        while (next >= 0 && !enabled(next)) next--;
        break;
      default:
        return;
    }
    e.preventDefault();
    if (next < 0) return;
    buttons[next].focus();
    options.pick(next);
  };
}
