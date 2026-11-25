/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { reportUnhandledError } from '@molstar/core/reactive/config';
import type { Observable } from '@molstar/core/reactive/observable';
import type { Child, Component } from './jsx-types.js';
import { runInScope, requireScope, toStream, ViewScope } from './scope.js';
import { bind, flush } from './watch.js';

const SVG_NS = 'http://www.w3.org/2000/svg';

/** SVG tags that are not also HTML tags; these are created in the SVG namespace. */
const SvgTags = new Set([
  'svg',
  'animate',
  'animateMotion',
  'animateTransform',
  'circle',
  'clipPath',
  'defs',
  'desc',
  'ellipse',
  'feBlend',
  'feColorMatrix',
  'feComponentTransfer',
  'feComposite',
  'feConvolveMatrix',
  'feDiffuseLighting',
  'feDisplacementMap',
  'feDistantLight',
  'feDropShadow',
  'feFlood',
  'feFuncA',
  'feFuncB',
  'feFuncG',
  'feFuncR',
  'feGaussianBlur',
  'feImage',
  'feMerge',
  'feMergeNode',
  'feMorphology',
  'feOffset',
  'fePointLight',
  'feSpecularLighting',
  'feSpotLight',
  'feTile',
  'feTurbulence',
  'filter',
  'foreignObject',
  'g',
  'image',
  'line',
  'linearGradient',
  'marker',
  'mask',
  'metadata',
  'mpath',
  'path',
  'pattern',
  'polygon',
  'polyline',
  'radialGradient',
  'rect',
  'set',
  'stop',
  'switch',
  'symbol',
  'text',
  'textPath',
  'tspan',
  'use',
  'view',
]); // prettier-ignore

/**
 * Tags that are both HTML and SVG (dom §4): SVG when made under an SVG parent (`svgContext`), when they hold SVG
 * content, or when placed as a static child of an SVG element (made again there).
 */
const SharedTags = new Set(['a', 'title', 'style', 'script']);

/**
 * Whether what renders now goes under an SVG parent (dom §4): set while a region, a `For` entry or `render` renders
 * into an element of the SVG namespace (not a `foreignObject`), kept while it renders into a fragment (its first
 * render, before its place is in the document).
 */
let svgContext = false;

/** Whether `host`'s children are SVG: an SVG element other than a `foreignObject`. */
function svgHost(host: Element): boolean {
  return host.namespaceURI === SVG_NS && host.localName !== 'foreignObject';
}

/**
 * Runs `fn` in `scope` with the namespace of `host`, where it renders (dom §4): an element sets it, anything else (a
 * fragment, nothing) keeps the current one.
 */
export function renderIn<T>(scope: ViewScope | undefined, host: Node | null | undefined, fn: () => T): T {
  const previous = svgContext;
  if (host instanceof Element) svgContext = svgHost(host);
  try {
    return runInScope(scope, fn);
  } finally {
    svgContext = previous;
  }
}

/** Whether a child holds an SVG element other than an `<svg>` root (an `<a>` around a `<text>` is SVG's own). */
function holdsSvg(child: unknown): boolean {
  if (child instanceof Element) return child.namespaceURI === SVG_NS && child.localName !== 'svg';
  if (child instanceof DocumentFragment) {
    for (let c = child.firstElementChild; c; c = c.nextElementSibling) if (holdsSvg(c)) return true;
    return false;
  }
  if (Array.isArray(child)) {
    for (const c of child) if (holdsSvg(c)) return true;
  }
  return false;
}

/** Set on a shared-tag HTML element whose props follow it (bindings, listeners, `ref`): made again, they would not. */
const Followed = Symbol('followed');
/** A shared-tag HTML element's static children (its `props.children`): made again with it, if it is. */
const StaticChildren = Symbol('staticChildren');

/**
 * A static child of an SVG element (its `props.children`: an element, or arrays of them) that is a shared-tag HTML
 * element, made again in the SVG namespace (dom §4): its attributes, delegated handlers and children moved, its own
 * static shared-tag children too. JSX makes children before their parent, so they could not know it. Only static
 * children: what a component or a control flow returns (a fragment) is referenced by it, and stays.
 */
function adoptStatic(child: unknown): unknown {
  if (Array.isArray(child)) {
    for (let i = 0; i < child.length; i++) child[i] = adoptStatic(child[i]);
    return child;
  }
  if (!(child instanceof Element) || child.namespaceURI === SVG_NS || !SharedTags.has(child.localName)) return child;
  const svg = document.createElementNS(SVG_NS, child.localName);
  for (const a of Array.from(child.attributes)) svg.setAttributeNS(a.namespaceURI, a.name, a.value);
  if ((child as any)[Handlers]) (svg as any)[Handlers] = (child as any)[Handlers];
  while (child.firstChild) svg.appendChild(child.firstChild);
  const kids = (child as any)[StaticChildren];
  if (kids !== undefined) adoptStaticIn(svg, kids);
  if ((child as any)[Followed]) {
    console.warn(
      `<${child.localName}> was made outside an SVG parent and is made again in one: its bound props, listeners and ref stay with the first element. Make it in a function child of the SVG element ({() => ...}), or bind on a child.`,
    );
  }
  return svg;
}

/** `adoptStatic` for the static children already in `parent` (an adopted element's). */
function adoptStaticIn(parent: Element, child: unknown): void {
  if (Array.isArray(child)) {
    for (const c of child) adoptStaticIn(parent, c);
    return;
  }
  const svg = adoptStatic(child);
  if (svg !== child && (child as Node).parentNode === parent) parent.replaceChild(svg as Node, child as Node);
}

/** Whether props hold something that follows the element made: a stream, an `on:` listener, a `ref`. */
function followsElement(props: Record<string, any>): boolean {
  for (const key in props) {
    if (key === 'children') continue;
    const value = props[key];
    if (key === 'ref' || key.startsWith('on:')) return true;
    // A delegated handler moves with the element; any other object (a stream, a style or classList record) binds.
    if (typeof value === 'function' && key.startsWith('on')) continue;
    if (!isStatic(value)) return true;
  }
  return false;
}

/** Properties that are read-only (or misbehave) on elements: always set as attributes. */
const AttributeOnly = new Set(['list', 'form', 'download', 'role', 'width', 'height']);

/**
 * Marks a component result that is data for its parent component rather than DOM (a `Match` for its `Switch`): `jsx`
 * passes it through as it is.
 */
export const Passthrough = Symbol('Passthrough');

/** Creates an element (or calls a component) with its props; children and bindings go to the current scope. */
export function createElement(type: string | Component<any>, props: Record<string, any> | null): Node {
  if (typeof type === 'function') {
    const result = type(props ?? {});
    if (result && typeof result === 'object' && Passthrough in result) return result as unknown as Node;
    return toNode(result);
  }
  const isSvg =
    SvgTags.has(type) || (SharedTags.has(type) && (svgContext || (props !== null && holdsSvg(props.children))));
  const el = isSvg ? document.createElementNS(SVG_NS, type) : document.createElement(type);
  if (props) {
    if (!isSvg && SharedTags.has(type)) {
      if (followsElement(props)) (el as any)[Followed] = true;
      if (props.children !== null && typeof props.children === 'object') (el as any)[StaticChildren] = props.children;
    }
    // A select's value only sticks once its options exist: apply it after the children.
    const deferred = type === 'select' && 'value' in props;
    applyProps(el, props, isSvg, type.includes('-'), deferred);
    let children = props.children;
    // Shared tags made before this parent (JSX makes children first) are made again in its namespace.
    if (isSvg && children !== null && typeof children === 'object' && svgHost(el)) children = adoptStatic(children);
    // Text alone: the element's text, without a node of our own to insert.
    if (typeof children === 'string' || typeof children === 'number') el.textContent = String(children);
    else if (children !== undefined) insert(el, children, null);
    if (deferred) applyProps(el, { value: props.value }, isSvg, false);
    if (typeof props.ref === 'function') props.ref(el);
  }
  return el;
}

/** Renders children into a fragment. */
export function Fragment(props: { children?: Child }): Node {
  const fragment = document.createDocumentFragment();
  insert(fragment, props.children, null);
  return fragment;
}

/** The child as one node: itself if it is one, otherwise a fragment holding it. */
export function toNode(child: Child): Node {
  if (child instanceof Node) return child;
  const fragment = document.createDocumentFragment();
  insert(fragment, child, null);
  return fragment;
}

/** Inserts a child before `before` (or at the end). Streams become regions that follow them. */
export function insert(parent: Node, child: Child, before: Node | null): void {
  if (child === null || child === undefined || typeof child === 'boolean') return;
  if (child instanceof Node) {
    parent.insertBefore(child, before);
  } else if (Array.isArray(child)) {
    for (const c of child) insert(parent, c, before);
  } else if (typeof child === 'function') {
    // Called in the parent's namespace (dom §4): a function child of an SVG element makes SVG.
    if (parent instanceof Element) {
      const previous = svgContext;
      svgContext = svgHost(parent);
      try {
        insert(parent, (child as () => Child)(), before);
      } finally {
        svgContext = previous;
      }
    } else insert(parent, (child as () => Child)(), before);
  } else if (typeof child === 'object') {
    const stream = toStream(child as Observable<Child>);
    if (stream) insertRegion(parent, stream, before);
    else reportUnhandledError(new Error(`Cannot render ${String(child)}`));
  } else {
    parent.insertBefore(document.createTextNode(String(child)), before);
  }
}

function isText(value: unknown): value is string | number | bigint | boolean | null | undefined {
  const t = typeof value;
  return value === null || t === 'undefined' || t === 'string' || t === 'number' || t === 'bigint' || t === 'boolean';
}

/**
 * A range between two comment anchors whose content follows a stream. Text values update one text node in place;
 * anything else replaces the range's content, rendered in a fresh child scope (the previous one is disposed).
 *
 * `memo`: a render function that is the one rendered last is not rendered again. The control-flow regions (`Show`,
 * `Switch`, `Dynamic`) emit the same function while their content stays, so that subscribing again after the scope was
 * paused and resumed (`ViewScope.setActive`) keeps what is on screen instead of rebuilding it.
 */
export function insertRegion(parent: Node, source: Observable<Child>, before: Node | null, memo = false): void {
  bind(source, regionWriter(parent, before, memo));
}

/**
 * The write of a region (`insertRegion`), for a caller that follows its own sources (`Show`, `Switch`, `Dynamic`): its
 * place is inserted now, before `before`; each call renders a value there.
 */
export function regionWriter(parent: Node, before: Node | null, memo = false): (value: Child) => void {
  const owner = requireScope('A reactive child');
  // Text first (the common case: a label, a count): one Text node is its own place. Anchors come only when the value
  // is nodes, around wherever the text was.
  let text: Text | undefined = document.createTextNode('');
  parent.insertBefore(text, before);
  let start: Comment | undefined;
  let end: Comment | undefined;
  let scope: ViewScope | undefined;
  let rendered: unknown;
  return (value) => {
    if (memo && typeof value === 'function' && value === rendered) return;
    rendered = value;
    if (isText(value)) {
      const data = value === null || value === undefined || typeof value === 'boolean' ? '' : String(value);
      if (text) {
        if (text.data !== data) text.data = data;
        return;
      }
      clearNodes();
      text = document.createTextNode(data);
      end!.parentNode!.insertBefore(text, end!);
      return;
    }
    if (!start) {
      start = document.createComment('');
      end = document.createComment('');
      const at = text!;
      at.parentNode!.insertBefore(start, at);
      at.parentNode!.insertBefore(end, at.nextSibling);
    }
    if (text) {
      text.remove();
      text = undefined;
    }
    clearNodes();
    const next = owner.child();
    scope = next;
    const fragment = document.createDocumentFragment();
    // In the namespace of the region's place (dom §4): content under an SVG element is SVG.
    renderIn(next, end!.parentNode, () => insert(fragment, value, null));
    end!.parentNode!.insertBefore(fragment, end!);
  };

  function clearNodes() {
    scope?.dispose();
    scope = undefined;
    if (start && end) removeBetween(start, end);
  }
}

/** Removes the nodes strictly between two siblings. */
export function removeBetween(start: Node, end: Node): void {
  const parent = start.parentNode;
  if (!parent) return;
  let node = start.nextSibling;
  while (node && node !== end) {
    const next = node.nextSibling;
    parent.removeChild(node);
    node = next;
  }
}

// --- Props ---------------------------------------------------------------------------------------

// --- Events --------------------------------------------------------------------------------------

/**
 * Bubbling events are delegated: an element keeps its handler (`onClick`), and each render root (a `render` container,
 * a `ViewElement`) listens once per type and walks from the target up to itself, calling the handlers it finds with
 * `currentTarget` as their own element; `stopPropagation` ends the walk as it ends propagation. A nested root goes
 * first and the next one (an outer root, then the document, for anything rendered outside a root) continues above it,
 * so nothing runs twice. The walk ends with a `flush`: what the handlers changed is rendered when the dispatch returns
 * (a discrete event, as a click dispatched by code, is on screen before the next line), not a microtask later. Events
 * that do not bubble (focus, blur, enter, leave) and those scrolling waits on (touch, wheel: they must stay passive)
 * stay native, and so does anything given as `on:type`. The pointer's move, over and out (pointer and mouse) stay native too: they come
 * at pointer rate over the whole page and a handful of elements take them, so a root walking each one's composed path
 * was all waste; a native listener runs only where one is.
 */
const DelegatedEvents = new Set([
  'click',
  'dblclick',
  'contextmenu',
  'input',
  'change',
  'beforeinput',
  'keydown',
  'keyup',
  'pointerdown',
  'pointerup',
  'pointercancel',
  'mousedown',
  'mouseup',
  'focusin',
  'focusout',
]);

const Handlers = Symbol('handlers');
const Walked = Symbol('walked');
const delegating = new WeakSet<Node>();
const eventTypes = new Map<string, string>();

/** Makes `root` handle the delegated events of what is rendered under it. */
export function delegateEvents(root: Node): void {
  if (delegating.has(root)) return;
  delegating.add(root);
  for (const type of DelegatedEvents) {
    root.addEventListener(type, dispatchDelegated);
    // An event made not to bubble (`new Event('input')` dispatched by code) never reaches the root's bubbling
    // listener: caught on its way down, it goes to the target's own handler only, as a native listener would.
    root.addEventListener(type, dispatchToTarget, true);
  }
}

function dispatchToTarget(e: Event) {
  if (e.bubbles || (e as any)[Walked]) return;
  (e as any)[Walked] = e.currentTarget;
  const target = e.target as Node;
  const handler = (target as any)[Handlers]?.[e.type];
  if (!handler) return;
  Object.defineProperty(e, 'currentTarget', { configurable: true, value: target });
  try {
    handler.call(target, e);
  } finally {
    delete (e as any).currentTarget;
    flush();
  }
}

function dispatchDelegated(e: Event) {
  if (!e.bubbles) return;
  const root = e.currentTarget as Node;
  // The path as it was when the event was dispatched: a handler that removes its own subtree (closing the Show it sits
  // in), or a detached inner root, does not cut the walk short, and the walk crosses shadow roots.
  const path = e.composedPath();
  const stop = path.indexOf(root);
  if (stop < 0) return;
  const from = (e as any)[Walked] as Node | undefined;
  // An inner root handled the path up to itself: continue above it.
  const at = from ? path.indexOf(from) : -1;
  let overridden = false;
  try {
    for (let i = from ? (at < 0 ? stop + 1 : at + 1) : 0; i <= stop; i++) {
      const node = path[i] as Node;
      const handler = (node as any)[Handlers]?.[e.type];
      if (!handler) continue;
      Object.defineProperty(e, 'currentTarget', { configurable: true, value: node });
      overridden = true;
      handler.call(node, e);
      if (e.cancelBubble) break;
    }
  } finally {
    (e as any)[Walked] = root;
    if (overridden) {
      delete (e as any).currentTarget;
      flush();
    }
  }
}

function setHandler(el: Element, type: string, handler: (e: Event) => void) {
  if (!DelegatedEvents.has(type)) {
    el.addEventListener(type, handler);
    return;
  }
  ((el as any)[Handlers] ??= {})[type] = handler;
  // Anything rendered outside a root is still covered (the document continues where the roots stopped).
  if (typeof document !== 'undefined') delegateEvents(document);
}

/** A value that is never a stream: written once, no binding (the common case, and the cheap one). */
const isStatic = (value: unknown) => value === null || (typeof value !== 'object' && typeof value !== 'function');

function applyProps(el: Element, props: Record<string, any>, isSvg: boolean, isCustom: boolean, skipValue = false) {
  let classes: ClassBinding | undefined;
  for (const key in props) {
    const value = props[key];
    if (key === 'children' || key === 'ref' || value === undefined || (skipValue && key === 'value')) continue;
    // A plain class with no flags: the attribute, without a binding.
    if (key === 'class' && typeof value === 'string' && props.classList === undefined) {
      if (value) el.setAttribute('class', value);
      continue;
    }
    if (key.startsWith('on:')) {
      el.addEventListener(key.slice(3), value);
    } else if (key.length > 2 && key.startsWith('on') && key.charCodeAt(2) < 91 && typeof value === 'function') {
      let type = eventTypes.get(key);
      if (type === undefined) eventTypes.set(key, (type = key.slice(2).toLowerCase()));
      setHandler(el, type, value);
    } else if (key === 'class' || key === 'className') {
      classes ??= new ClassBinding(el);
      bind(value, (v: string | undefined) => classes!.setBase(v));
    } else if (key === 'classList') {
      classes ??= new ClassBinding(el);
      for (const name of Object.keys(value)) {
        bind(value[name], (on: boolean | undefined) => classes!.setFlag(name, !!on));
      }
    } else if (key === 'style') {
      applyStyle(el as HTMLElement, value);
    } else if (key.startsWith('attr:')) {
      const name = key.slice(5);
      bind(value, (v) => setAttribute(el, name, v));
    } else if (key.startsWith('prop:')) {
      const name = key.slice(5);
      bind(value, (v) => {
        (el as any)[name] = v;
      });
    } else if (isSvg || key.includes('-')) {
      if (isStatic(value)) setAttribute(el, key, value);
      else bind(value, (v) => setAttribute(el, key, v));
    } else if (isCustom) {
      // Custom elements take models, lists and callbacks as properties.
      bind(value, (v) => {
        (el as any)[key] = v;
      });
    } else if (key in el && !AttributeOnly.has(key)) {
      if (isStatic(value)) {
        (el as any)[key] = value === null ? '' : value;
        if (value === null || value === false) el.removeAttribute(key);
        continue;
      }
      bind(value, (v) => {
        (el as any)[key] = v === null || v === undefined ? '' : v;
        if (v === null || v === undefined || v === false) el.removeAttribute(key);
      });
    } else if (isStatic(value)) {
      setAttribute(el, key, value);
    } else {
      bind(value, (v) => setAttribute(el, key, v));
    }
  }
}

function setAttribute(el: Element, name: string, value: unknown) {
  // ARIA states are the strings "true"/"false", not present/absent.
  if (typeof value === 'boolean' && name.startsWith('aria-')) el.setAttribute(name, String(value));
  else if (value === null || value === undefined || value === false) el.removeAttribute(name);
  else el.setAttribute(name, value === true ? '' : String(value));
}

/** `class` and `classList` together: the class attribute is the base classes plus the flags that are on. */
class ClassBinding {
  #base = '';
  readonly #flags = new Map<string, boolean>();
  constructor(readonly el: Element) {}

  setBase(value: string | undefined) {
    this.#base = value ?? '';
    this.#sync();
  }

  setFlag(name: string, on: boolean) {
    this.#flags.set(name, on);
    this.#sync();
  }

  #sync() {
    let value = this.#base;
    for (const [name, on] of this.#flags) if (on) value = value ? `${value} ${name}` : name;
    if (value) this.el.setAttribute('class', value);
    else this.el.removeAttribute('class');
  }
}

function applyStyle(el: HTMLElement | SVGElement, value: unknown) {
  if (typeof value === 'string' || toStream(value as string)) {
    bind(value as string | Observable<string | undefined>, (v) => {
      if (v) el.setAttribute('style', v);
      else el.removeAttribute('style');
    });
    return;
  }
  const styles = value as Record<string, unknown>;
  for (const name of Object.keys(styles)) {
    // Custom properties and kebab-case names go through setProperty; camelCase through the style object.
    const property = name.startsWith('--') || name.includes('-') ? name : toKebab(name);
    const px = !property.startsWith('--') && !Unitless.has(property);
    bind(styles[name], (v) => {
      if (v === null || v === undefined || v === '') el.style.removeProperty(property);
      else el.style.setProperty(property, typeof v === 'number' && px && v !== 0 ? `${v}px` : String(v));
    });
  }
}

/** Properties whose numbers have no unit; other numbers are pixels. */
const Unitless = new Set([
  'opacity',
  'z-index',
  'flex',
  'flex-grow',
  'flex-shrink',
  'order',
  'line-height',
  'font-weight',
  'zoom',
  'aspect-ratio',
  'scale',
  'fill-opacity',
  'stroke-opacity',
  'stroke-width',
  'grid-row',
  'grid-column',
  'tab-size',
  'orphans',
  'widows',
  'columns',
  'column-count',
]); // prettier-ignore

function toKebab(name: string) {
  return name.replace(/[A-Z]/g, (c) => `-${c.toLowerCase()}`);
}

// --- Roots ---------------------------------------------------------------------------------------

/**
 * Renders into a container, after what is already there, in a new root scope (or, with `parent`, a child of it: it
 * pauses and resumes with it, `ViewScope.setActive`, which is how a view mounted in a hidden panel stays quiet). Returns a
 * function that disposes the scope and removes what was rendered.
 */
export function render(fn: () => Child, container: Node, parent?: ViewScope): () => void {
  delegateEvents(container);
  const scope = parent ? parent.child() : new ViewScope();
  const start = document.createComment('');
  const end = document.createComment('');
  const fragment = document.createDocumentFragment();
  fragment.append(start);
  renderIn(scope, container, () => insert(fragment, fn(), null));
  fragment.append(end);
  container.appendChild(fragment);
  return () => {
    scope.dispose();
    removeBetween(start, end);
    start.remove();
    end.remove();
  };
}
