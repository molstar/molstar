/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { BehaviorSubject, Observable, ReactiveModel, Scope, State, Subject, Subscription } from '../index.js';

// The model lifecycle and its ownership (reactive §7).

describe('Scope', () => {
  it('ends what it owns last added first, attempts all, rethrows the first failure', () => {
    const order: string[] = [];
    const scope = new Scope();
    scope.own(() => order.push('a'));
    scope.own({ dispose: () => order.push('b') });
    scope.own({
      unsubscribe: () => {
        order.push('c');
        throw new Error('c failed');
      },
    });
    scope.own(() => {
      order.push('d');
      throw new Error('d failed');
    });
    expect(() => scope.dispose()).toThrow('d failed');
    expect(order).toEqual(['d', 'c', 'b', 'a']);
    scope.own(() => order.push('late'));
    expect(order.at(-1)).toBe('late');
  });

  it('a subscription that ends first leaves it', () => {
    const scope = new Scope();
    const s = new Subject<number>();
    for (let i = 0; i < 50; i++) {
      const sub = s.subscribe(() => {});
      scope.add(sub);
      sub.unsubscribe();
    }
    expect((scope as any)._finalizers).toBeUndefined();
  });
});

describe('ReactiveModel', () => {
  it('gives each mount an abort signal that ends with it', () => {
    class Example extends ReactiveModel {
      read = () => this.mountSignal;
      inMount?: AbortSignal;
      protected override onMount() {
        this.inMount = this.mountSignal;
      }
    }
    const model = new Example();
    expect(model.read().aborted).toBe(true);
    model.mount();
    const first = model.inMount!;
    expect(first.aborted).toBe(false);
    expect(model.read()).toBe(first);
    model.unmount();
    expect(first.aborted).toBe(true);
    model.mount();
    const second = model.read();
    expect(second).not.toBe(first);
    model.retire();
    expect(second.aborted).toBe(true);
  });

  it('mounts once per cycle, ends what the mount registered in reverse, and remounts', () => {
    const order: string[] = [];
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    class Example extends ReactiveModel {
      protected override onMount(): void {
        order.push('mount');
        this.registerDispose(() => order.push('first'));
        this.registerDispose(() => order.push('second'));
      }
    }
    const model = new Example();
    model.mount();
    model.mount();
    model.mount();
    model.unmount();
    expect(order).toEqual(['mount', 'second', 'first']);
    expect(warn).toHaveBeenCalledTimes(1);
    model.mount();
    model.unmount();
    expect(order).toEqual(['mount', 'second', 'first', 'mount', 'second', 'first']);
    expect(warn).toHaveBeenCalledTimes(1);
  });

  it('subscribe: the current value, then each change, until unmount', () => {
    const source = new State(1);
    const seen: number[] = [];
    class Example extends ReactiveModel {
      protected override onMount(): void {
        this.subscribe(source, (v) => seen.push(v));
      }
    }
    const model = new Example();
    model.mount();
    source.set(2);
    model.unmount();
    source.set(3);
    expect(seen).toEqual([1, 2]);
    expect(source.observed).toBe(false);
  });

  it('listen: each change, no current value, until unmount; one object per listener on a state', () => {
    const source = new State(1);
    const subject = new BehaviorSubject('a');
    const seen: unknown[] = [];
    class Example extends ReactiveModel {
      protected override onMount(): void {
        this.listen(source, (v) => seen.push(v));
        this.listen(source.asReadonly, (v) => seen.push(`view ${v}`));
        this.listen(subject, (v) => seen.push(v));
      }
    }
    const model = new Example();
    model.mount();
    source.set(2);
    subject.next('b');
    model.unmount();
    source.set(3);
    subject.next('c');
    expect(seen).toEqual([2, 'view 2', 'b']);
    expect(source.observed || subject.observed).toBe(false);
  });

  it('own and ownPermanently end with the mount and with retire', () => {
    const order: string[] = [];
    class Example extends ReactiveModel {
      constructor() {
        super();
        this.ownPermanently({ dispose: () => order.push('permanent') });
      }
      protected override onMount(): void {
        this.own(new Subscription(() => order.push('subscription')));
        this.own({ dispose: () => order.push('disposable') });
      }
    }
    const model = new Example();
    model.mount();
    model.unmount();
    expect(order).toEqual(['disposable', 'subscription']);
    model.retire();
    expect(order).toEqual(['disposable', 'subscription', 'permanent']);
  });

  it('a synchronous value that unmounts the model ends the subscription', () => {
    let active = 0;
    const source = new Observable<number>((s) => {
      active++;
      s.next(1);
      return () => active--;
    });
    class Example extends ReactiveModel {
      protected override onMount(): void {
        this.subscribe(source, () => this.unmount());
      }
    }
    const model = new Example();
    model.mount();
    expect(active).toBe(0);
    model.mount();
    expect(active).toBe(0);
  });

  it('a synchronous value that retires the model ends a permanent subscription', () => {
    let active = 0;
    const source = new Observable<number>((s) => {
      active++;
      s.next(1);
      return () => active--;
    });
    class Example extends ReactiveModel {
      subscribeNow(): void {
        this.subscribePermanently(source, () => this.retire());
      }
    }
    const model = new Example();
    model.subscribeNow();
    expect(active).toBe(0);
    expect(() => model.subscribeNow()).toThrow('Cannot subscribe permanently after retirement');
  });

  it('an unmount during onMount leaves nothing live', () => {
    const closed = new BehaviorSubject(true);
    const other = new BehaviorSubject(0);
    class Example extends ReactiveModel {
      protected override onMount(): void {
        this.subscribe(closed, (v) => {
          if (v) this.unmount();
        });
        this.subscribe(other, () => {});
      }
    }
    new Example().mount();
    expect(closed.observed || other.observed).toBe(false);
  });

  it('registers, subscribes, listens and owns only in onMount', () => {
    class Example extends ReactiveModel {
      register = () => this.registerDispose(() => {});
      sub = () => this.subscribe(new Subject(), () => {});
      lis = () => this.listen(new Subject(), () => {});
      ownIt = () => this.own(() => {});
    }
    const model = new Example();
    model.mount();
    for (const f of [model.register, model.sub, model.lis, model.ownIt]) expect(f).toThrow(/outside onMount/);
  });

  it('rolls back a failed mount, attempting every cleanup, and allows a later mount', () => {
    const order: string[] = [];
    class Example extends ReactiveModel {
      attempts = 0;
      protected override onMount(): void {
        this.attempts++;
        this.registerDispose(() => order.push('first'));
        this.registerDispose(() => {
          order.push('second');
          if (this.attempts === 1) throw new Error('cleanup failed');
        });
        if (this.attempts === 1) throw new Error('mount failed');
        order.push('mounted');
      }
    }
    const model = new Example();
    expect(() => model.mount()).toThrow('mount failed');
    expect(order).toEqual(['second', 'first']);
    model.mount();
    model.unmount();
    expect(order).toEqual(['second', 'first', 'mounted', 'second', 'first']);
  });

  it('runs the cleanup when onUnmount fails, and the permanent cleanup when retiring fails', () => {
    const order: string[] = [];
    class Example extends ReactiveModel {
      constructor() {
        super();
        this.registerPermanentDispose(() => order.push('permanent'));
      }
      protected override onMount(): void {
        this.registerDispose(() => order.push('mount cleanup'));
      }
      protected override onUnmount(): void {
        order.push('unmount');
        throw new Error('unmount failed');
      }
      protected override onRetire(): void {
        order.push('retire');
        throw new Error('retire failed');
      }
    }
    const model = new Example();
    model.mount();
    expect(() => model.retire()).toThrow('unmount failed');
    expect(order).toEqual(['unmount', 'mount cleanup', 'retire', 'permanent']);
    expect(() => model.retire()).not.toThrow();
    model.mount();
    expect(order.length).toBe(4);
  });

  it('a remover runs its dispose once, now; the scope does not run it again', () => {
    const order: string[] = [];
    class Example extends ReactiveModel {
      removeChild: (() => void) | undefined;
      protected override onMount(): void {
        this.registerDispose(() => order.push('retained'));
        this.removeChild = this.registerDispose(() => order.push('child'));
      }
    }
    const model = new Example();
    model.mount();
    model.removeChild!();
    model.unmount();
    model.removeChild!();
    expect(order).toEqual(['child', 'retained']);
  });

  it('keeps the lifecycle methods with the base class', () => {
    class Overrides extends ReactiveModel {
      override mount(): void {}
    }
    expect(() => new Overrides()).toThrow(/cannot be overridden/);
    class OverridesRetire extends ReactiveModel {
      override retire(): void {}
    }
    expect(() => new OverridesRetire()).toThrow(/cannot be overridden/);
    const model = new (class extends ReactiveModel {})();
    expect(() => {
      (model as any).mount = () => {};
    }).toThrow();
  });
});
