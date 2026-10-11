/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Observable } from 'rxjs';
import type { MaybeReactive } from './scope.js';

/**
 * What can be rendered: DOM nodes, text, nothing, lists, a stream or state of any of these (the DOM follows it), or a
 * function rendered in its own scope.
 */
export type Child =
  | Node
  | string
  | number
  | bigint
  | boolean
  | null
  | undefined
  | Observable<Child>
  | { readonly asReadonly: Observable<Child> }
  | readonly Child[]
  | (() => Child);

export type Component<P = {}> = (props: P) => Child;

/** Event props in the usual camel case; the runtime lowercases the name after `on`. */
interface EventNames {
  onClick: 'click';
  onDblClick: 'dblclick';
  onContextMenu: 'contextmenu';
  onAuxClick: 'auxclick';
  onKeyDown: 'keydown';
  onKeyUp: 'keyup';
  onInput: 'input';
  onChange: 'change';
  onBeforeInput: 'beforeinput';
  onFocus: 'focus';
  onBlur: 'blur';
  onFocusIn: 'focusin';
  onFocusOut: 'focusout';
  onPointerDown: 'pointerdown';
  onPointerUp: 'pointerup';
  onPointerMove: 'pointermove';
  onPointerEnter: 'pointerenter';
  onPointerLeave: 'pointerleave';
  onPointerCancel: 'pointercancel';
  onMouseDown: 'mousedown';
  onMouseUp: 'mouseup';
  onMouseMove: 'mousemove';
  onMouseEnter: 'mouseenter';
  onMouseLeave: 'mouseleave';
  onWheel: 'wheel';
  onScroll: 'scroll';
  onSubmit: 'submit';
  onDragStart: 'dragstart';
  onDragEnter: 'dragenter';
  onDragOver: 'dragover';
  onDragLeave: 'dragleave';
  onDrop: 'drop';
  onDragEnd: 'dragend';
  onCopy: 'copy';
  onCut: 'cut';
  onPaste: 'paste';
  onToggle: 'toggle';
  onAnimationEnd: 'animationend';
  onTransitionEnd: 'transitionend';
}

export type EventProps<E> = {
  [K in keyof EventNames]?: (e: HTMLElementEventMap[EventNames[K]] & { currentTarget: E }) => void;
};

/** Props every element takes. */
export interface CommonProps<E extends Element> {
  children?: Child;
  /** Receives the element once its props are set. */
  ref?: (el: E) => void;
  class?: MaybeReactive<string | undefined>;
  /** Classes toggled by flags; combined with `class`. */
  classList?: Readonly<Record<string, MaybeReactive<boolean | undefined>>>;
  style?: MaybeReactive<string | undefined> | Readonly<Record<string, MaybeReactive<string | number | undefined>>>;
  /** Any event by its exact name: `on:my-event`. */
  [event: `on:${string}`]: ((e: Event) => void) | undefined;
  /** Forces an attribute: `attr:popovertarget`. */
  [attr: `attr:${string}`]: MaybeReactive<string | number | boolean | undefined>;
  /** Forces a property: `prop:value`. */
  [prop: `prop:${string}`]: unknown;
}

/**
 * The element's own properties (not methods or `on…` handler properties), each optionally reactive; `null` or
 * `undefined` clears the attribute.
 */
type PropertyKeys<E> = {
  [K in keyof E]: E[K] extends (...args: any[]) => any ? never : K extends `on${string}` ? never : K;
}[keyof E];

export type ElementProps<E extends Element> = CommonProps<E> &
  EventProps<E> & {
    [K in Exclude<PropertyKeys<E>, keyof CommonProps<E> | 'className' | 'classList' | 'style'>]?: MaybeReactive<
      E[K] | null | undefined
    >;
  };

/** SVG is typed loosely: attributes by their SVG names (`viewBox`, `stroke-width`). */
export type SVGProps = CommonProps<SVGElement> & EventProps<SVGElement> & { [attr: string]: unknown };

type HTMLElements = { [K in keyof HTMLElementTagNameMap]: ElementProps<HTMLElementTagNameMap[K]> };
type SVGElements = { [K in Exclude<keyof SVGElementTagNameMap, keyof HTMLElementTagNameMap>]: SVGProps };

/**
 * Custom elements (`msp-*`) add themselves here by module augmentation:
 * `declare module '../dom/jsx-types.js' { interface MspElements { 'msp-x': MspElementProps<MspX> } }`.
 */
export interface MspElements {}

export declare namespace JSX {
  type Element = Node;
  type ElementType = keyof IntrinsicElements | ((props: any) => Child);
  interface ElementChildrenAttribute {
    children: {};
  }
  interface IntrinsicAttributes {}
  interface IntrinsicElements extends HTMLElements, SVGElements, MspElements {}
}
