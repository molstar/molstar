/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ReactiveModel } from '@molstar/core/reactive/model';
import type { Child, CommonProps } from './jsx-types.js';
import { DomHost } from './host.js';
import { delegateEvents, insert } from './render.js';
import { runInScope, ViewScope } from './scope.js';

/** `HTMLElement` where it exists, so modules defining elements can be imported outside a browser (Node, workers). */
const ElementBase = (typeof HTMLElement === 'undefined' ? class {} : HTMLElement) as typeof HTMLElement;

/** Where an element is: its document, or the shadow root it is in. */
export type ElementRoot = Document | ShadowRoot;

/**
 * A custom element that renders a view of its `model` into its own light DOM, once, with bindings doing the updates
 * (dom §7).
 *
 * - Rendering starts when the element is connected and has a model; setting another model renders again.
 * - Disconnecting disposes the rendering on the next microtask, unless the element is connected again before then:
 *   moving an element (between containers, into a dialog) keeps its DOM and bindings.
 * - `onConnected(root)` runs at each connect with where it is (a document or a shadow root; the application adopts its
 *   styles there), and `onReleased(root)` once it has left that root and not come back by the next microtask.
 * - The element does not own its model: whoever created the model mounts and retires it. An element that creates its
 *   own model (`createModel`) mounts it while rendered and unmounts it when the rendering is disposed.
 */
export abstract class ViewElement<M = unknown> extends ElementBase {
  private _model: M | undefined;
  private _owned = false;
  private _scope: ViewScope | undefined;
  private _pendingDispose = false;
  private _root: ElementRoot | undefined;

  constructor() {
    super();
    // A `model` set before the element was defined is an own property that hides the accessor: take it over.
    if (Object.hasOwn(this, 'model')) {
      const model = (this as any).model as M;
      delete (this as any).model;
      this.model = model;
    }
  }

  get model(): M | undefined {
    return this._model;
  }

  set model(model: M | undefined) {
    if (model === this._model) return;
    this._unrender();
    this._model = model;
    this._owned = false;
    if (this.isConnected) this._render();
  }

  /** Renders the view of the model; called once per model while connected. */
  protected abstract render(model: M): Child;

  /** Creates a model when none was given; the element then owns it. */
  protected createModel?(): M;

  /** Connected in `root` (each connect): the application's hook (styles, a root class). */
  protected onConnected?(root: ElementRoot): void;

  /** Left `root` for good (not connected there again by the next microtask). */
  protected onReleased?(root: ElementRoot): void;

  connectedCallback(): void {
    const node = this.getRootNode();
    const root = typeof ShadowRoot !== 'undefined' && node instanceof ShadowRoot ? node : this.ownerDocument;
    if (root !== this._root) {
      if (this._root) this.onReleased?.(this._root);
      this._root = root;
    }
    this.onConnected?.(root);
    if (this._pendingDispose) {
      this._pendingDispose = false;
      return;
    }
    this._render();
  }

  disconnectedCallback(): void {
    DomHost.queueMicrotask(() => {
      if (this.isConnected || !this._root) return;
      const root = this._root;
      this._root = undefined;
      this.onReleased?.(root);
    });
    if (!this._scope) return;
    this._pendingDispose = true;
    DomHost.queueMicrotask(() => {
      if (!this._pendingDispose) return;
      this._pendingDispose = false;
      this._unrender();
    });
  }

  private _render(): void {
    if (this._scope) return;
    if (this._model === undefined && this.createModel) {
      this._model = this.createModel();
      this._owned = true;
    }
    const model = this._model;
    if (model === undefined) return;
    const scope = new ViewScope();
    this._scope = scope;
    if (this._owned && model instanceof ReactiveModel) {
      model.mount();
      scope.onCleanup(() => model.unmount());
    }
    delegateEvents(this);
    const content = document.createDocumentFragment();
    runInScope(scope, () => insert(content, this.render(model), null));
    this.replaceChildren(content);
  }

  private _unrender(): void {
    this._pendingDispose = false;
    const scope = this._scope;
    if (!scope) return;
    this._scope = undefined;
    scope.dispose();
    this.replaceChildren();
  }
}

/** JSX props of a `ViewElement`: the common element props plus its model. */
export type ViewElementProps<E extends ViewElement<any>> = CommonProps<E> & { model?: E['model'] };

export interface ViewElementClass {
  readonly tag: string;
  new (): ViewElement<any>;
}

/** Registers elements under their tags; an already registered tag is left as it is. */
export function defineElements(elements: readonly ViewElementClass[]): void {
  for (const element of elements) {
    if (!customElements.get(element.tag)) customElements.define(element.tag, element);
  }
}
