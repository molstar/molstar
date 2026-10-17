/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * The writes of a `watch` render (rendering.md §2): each compares with what the DOM has and writes only when it
 * differs, so a render that runs for any change of its sources touches only what changed. Nothing is kept or allocated
 * per write, and nothing read forces layout (text, attributes, classes, inline styles).
 */

/** A text node's text. */
export function setText(node: Text, value: string): void {
  if (node.data !== value) node.data = value;
}

/**
 * An attribute: a string, a number, true (present, empty), false/null/undefined (removed). ARIA states are the strings
 * "true" and "false": give them as strings.
 */
export function setAttr(el: Element, name: string, value: string | number | boolean | null | undefined): void {
  if (value === undefined || value === null || value === false) {
    if (el.hasAttribute(name)) el.removeAttribute(name);
    return;
  }
  const v = value === true ? '' : String(value);
  if (el.getAttribute(name) !== v) el.setAttribute(name, v);
}

/** A class, on or off (`classList.toggle` with a force writes nothing when it is already so). */
export function setClass(el: Element, name: string, on: boolean): void {
  el.classList.toggle(name, on);
}

/** An inline style property or custom property ('', null or undefined: removed). */
export function setStyle(el: HTMLElement | SVGElement, name: string, value: string | null | undefined): void {
  const style = el.style;
  if (!value) {
    if (style.getPropertyValue(name)) style.removeProperty(name);
    return;
  }
  if (style.getPropertyValue(name) !== value) style.setProperty(name, value);
}

/** Shown, or `display: none` (the attribute `hidden` would lose to a class's `display`). */
export function setShown(el: HTMLElement | SVGElement, on: boolean): void {
  const display = on ? '' : 'none';
  if (el.style.display !== display) el.style.display = display;
}
