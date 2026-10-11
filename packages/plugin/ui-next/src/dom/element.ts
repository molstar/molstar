/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ReactiveModel } from '../reactive/reactive-model.js';
import type { Child, CommonProps } from './jsx-types.js';
import { insert } from './render.js';
import { runInScope, Scope } from './scope.js';

/** `HTMLElement` where it exists, so modules defining elements can be imported outside a browser (Node, workers). */
const ElementBase = (typeof HTMLElement === 'undefined' ? class {} : HTMLElement) as typeof HTMLElement;

/**
 * A custom element that renders a view of its `model` into its own light DOM, once, with bindings doing the updates.
 *
 * - Rendering starts when the element is connected and has a model; setting another model renders again.
 * - Disconnecting disposes the rendering on the next microtask, unless the element is connected again before then:
 *   moving an element (between docks, into a dialog) keeps its DOM and bindings.
 * - The element does not own its model: whoever created the model mounts and retires it. An element that creates its
 *   own model (`createModel`) mounts it while rendered and unmounts it when the rendering is disposed.
 */
export abstract class MspElement<M = unknown> extends ElementBase {
  #model: M | undefined;
  #owned = false;
  #scope: Scope | undefined;
  #pendingDispose = false;

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
    return this.#model;
  }

  set model(model: M | undefined) {
    if (model === this.#model) return;
    this.#unrender();
    this.#model = model;
    this.#owned = false;
    if (this.isConnected) this.#render();
  }

  /** Renders the view of the model; called once per model while connected. */
  protected abstract render(model: M): Child;

  /** Creates a model when none was given; the element then owns it. */
  protected createModel?(): M;

  connectedCallback(): void {
    if (this.#pendingDispose) {
      this.#pendingDispose = false;
      return;
    }
    this.#render();
  }

  disconnectedCallback(): void {
    if (!this.#scope) return;
    this.#pendingDispose = true;
    queueMicrotask(() => {
      if (!this.#pendingDispose) return;
      this.#pendingDispose = false;
      this.#unrender();
    });
  }

  #render(): void {
    if (this.#scope) return;
    if (this.#model === undefined && this.createModel) {
      this.#model = this.createModel();
      this.#owned = true;
    }
    const model = this.#model;
    if (model === undefined) return;
    const scope = new Scope();
    this.#scope = scope;
    if (this.#owned && model instanceof ReactiveModel) {
      model.mount();
      scope.onCleanup(() => model.unmount());
    }
    const content = document.createDocumentFragment();
    runInScope(scope, () => insert(content, this.render(model), null));
    this.replaceChildren(content);
  }

  #unrender(): void {
    this.#pendingDispose = false;
    const scope = this.#scope;
    if (!scope) return;
    this.#scope = undefined;
    scope.dispose();
    this.replaceChildren();
  }
}

/** JSX props of an `MspElement`: the common element props plus its model. */
export type MspElementProps<E extends MspElement<any>> = CommonProps<E> & { model?: E['model'] };

export interface MspElementClass {
  readonly tag: string;
  new (): MspElement<any>;
}

/** Registers elements under their tags; an already registered tag is left as it is. */
export function defineElements(elements: readonly MspElementClass[]): void {
  for (const element of elements) {
    if (!customElements.get(element.tag)) customElements.define(element.tag, element);
  }
}
