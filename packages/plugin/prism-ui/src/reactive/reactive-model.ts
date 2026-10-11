/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Observable, Subscription } from 'rxjs';

import { createDisposableGroup, type DisposableGroup } from './disposable-group.js';

export abstract class ReactiveModel {
  #permanentGroup = createDisposableGroup();
  #mountGroup: DisposableGroup | undefined;
  #mounting = false;
  #mounted = false;
  #retired = false;
  // An unmount asked for while onMount runs, carried out when it returns.
  #unmountAfterMount = false;
  // A second mount while mounted means two owners; said once.
  #warnedTwoOwners = false;

  constructor() {
    if (
      this.mount !== ReactiveModel.prototype.mount ||
      this.unmount !== ReactiveModel.prototype.unmount ||
      this.retire !== ReactiveModel.prototype.retire
    ) {
      throw new Error('ReactiveModel lifecycle methods cannot be overridden; use the protected lifecycle hooks');
    }
    Object.defineProperties(this, {
      mount: { configurable: false, enumerable: false, value: () => this.#mount(), writable: false },
      unmount: { configurable: false, enumerable: false, value: () => this.#unmount(), writable: false },
      retire: { configurable: false, enumerable: false, value: () => this.#retire(), writable: false },
    });
  }

  mount(): void {
    this.#mount();
  }

  #mount(): void {
    if (this.#retired) {
      return;
    }
    if (this.#mounted) {
      if (!this.#warnedTwoOwners) {
        this.#warnedTwoOwners = true;
        console.warn(
          `${this.constructor.name}: mounted again while mounted. A model has one owner: the first unmount will tear it down under the other.`,
        );
      }
      return;
    }
    const mountGroup = createDisposableGroup();
    this.#mounted = true;
    this.#mountGroup = mountGroup;
    this.#mounting = true;
    try {
      this.onMount();
    } catch (error) {
      this.#mounted = false;
      if (this.#mountGroup === mountGroup) {
        this.#mountGroup = undefined;
      }
      try {
        mountGroup.dispose();
      } catch {
        // The mount hook failed first; cleanup remains best effort and must not hide that failure.
      }
      throw error;
    } finally {
      this.#mounting = false;
    }
    if (this.#unmountAfterMount) {
      this.#unmountAfterMount = false;
      this.#unmount();
    }
  }

  unmount(): void {
    this.#unmount();
  }

  #unmount(): void {
    if (!this.#mounted) {
      return;
    }
    if (this.#mounting) {
      this.#unmountAfterMount = true;
      return;
    }
    this.#mounted = false;
    const mountGroup = this.#mountGroup;
    this.#mountGroup = undefined;
    let failed = false;
    let failure: unknown;
    try {
      this.onUnmount();
    } catch (error) {
      failed = true;
      failure = error;
    }
    try {
      mountGroup?.dispose();
    } catch (error) {
      if (!failed) {
        failed = true;
        failure = error;
      }
    }
    if (failed) {
      throw failure;
    }
  }

  /** Irreversibly ends this model. Only its final owner calls this; React cleanup uses unmount. */
  retire(): void {
    this.#retire();
  }

  #retire(): void {
    if (this.#retired) {
      return;
    }
    this.#retired = true;
    let failed = false;
    let failure: unknown;
    try {
      this.#unmount();
    } catch (error) {
      failed = true;
      failure = error;
    }
    try {
      this.onRetire();
    } catch (error) {
      if (!failed) {
        failed = true;
        failure = error;
      }
    }
    try {
      this.#permanentGroup.dispose();
    } catch (error) {
      if (!failed) {
        failed = true;
        failure = error;
      }
    }
    if (failed) {
      throw failure;
    }
  }

  protected registerDispose(dispose: () => void): () => void {
    if (!this.#mounting || this.#mountGroup === undefined) {
      throw new Error('Cannot register dispose outside onMount');
    }
    return this.#mountGroup.registerDispose(dispose);
  }

  protected registerPermanentDispose(dispose: () => void): () => void {
    if (this.#retired) {
      throw new Error('Cannot register permanent dispose after retirement');
    }
    return this.#permanentGroup.registerDispose(dispose);
  }

  protected subscribe<T>(source: Observable<T>, action: (value: T) => void): void {
    this.#subscribeOwned((dispose) => this.registerDispose(dispose), source, action);
  }

  protected subscribePermanently<T>(source: Observable<T>, action: (value: T) => void): void {
    if (this.#retired) {
      throw new Error('Cannot subscribe permanently after retirement');
    }
    this.#subscribeOwned((dispose) => this.registerPermanentDispose(dispose), source, action);
  }

  #subscribeOwned<T>(
    register: (dispose: () => void) => () => void,
    source: Observable<T>,
    action: (value: T) => void,
  ): void {
    let subscription: Subscription | undefined;
    let released = false;
    const release = register(() => {
      released = true;
      subscription?.unsubscribe();
    });
    try {
      subscription = source.subscribe(action);
    } catch (error) {
      release();
      throw error;
    }
    if (released) {
      subscription.unsubscribe();
    }
  }

  protected onMount(): void {}

  protected onUnmount(): void {}

  protected onRetire(): void {}
}
