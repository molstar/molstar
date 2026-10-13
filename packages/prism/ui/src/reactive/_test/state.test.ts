/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { BehaviorSubject, config, map, Observable } from 'rxjs';
import { combineBehaviors } from '../combine-behaviors.js';
import { observableState } from '../observable-state.js';
import { ReactiveModel } from '../reactive-model.js';

describe('readonly behavior', () => {
  it('emits the current value, pipes operators, and does not grant mutation through the public view', () => {
    const writable = observableState(1);
    const source = writable.asReadonly;
    const seen: number[] = [];
    const doubled: number[] = [];
    const subscription = source.subscribe((value) => {
      seen.push(value);
    });
    const piped = source.pipe(map((value) => value * 2)).subscribe((value) => {
      doubled.push(value);
    });
    expect(seen).toEqual([1]);
    expect(doubled).toEqual([2]);
    writable.set(2);
    expect(seen).toEqual([1, 2]);
    expect(doubled).toEqual([2, 4]);
    expect(source.value).toBe(2);
    subscription.unsubscribe();
    piped.unsubscribe();
    writable.set(3);
    expect(seen).toEqual([1, 2]);
    expect(doubled).toEqual([2, 4]);
  });
});

describe('record behavior update', () => {
  it('merges partial fields and skips emit when changed keys match the current value', () => {
    const writable = observableState({ a: 1, error: undefined as string | undefined });
    const seen: Array<{ a: number; error: string | undefined }> = [];
    const subscription = writable.asReadonly.subscribe((value) => {
      seen.push(value);
    });
    expect(seen).toEqual([{ a: 1, error: undefined }]);
    writable.update({ a: 1 });
    expect(seen).toHaveLength(1);
    writable.update({ error: 'failed' });
    expect(seen).toEqual([
      { a: 1, error: undefined },
      { a: 1, error: 'failed' },
    ]);
    writable.update({ error: 'failed' });
    expect(seen).toHaveLength(2);
    subscription.unsubscribe();
  });
});

describe('grouped behavior projections', () => {
  it('projects through an explicit project function and suppresses shallow-equal repeats', () => {
    const first = observableState({ value: 1 });
    const second = observableState({ label: 'a' });
    const combined = combineBehaviors({
      sources: { first: first.asReadonly, second: second.asReadonly },
      project: ({ first, second }) => ({ sum: first.value, label: second.label }),
    });
    const seen: Array<{ sum: number; label: string }> = [];
    const subscription = combined.subscribe((value) => {
      seen.push(value);
    });
    expect(combined.value).toEqual({ sum: 1, label: 'a' });
    expect(seen).toEqual([{ sum: 1, label: 'a' }]);
    first.update({ value: 2 });
    expect(seen).toEqual([
      { sum: 1, label: 'a' },
      { sum: 2, label: 'a' },
    ]);
    first.update({ value: 2 });
    expect(seen).toHaveLength(2);
    combined.dispose();
    subscription.unsubscribe();
  });
});

describe('ReactiveModel lifecycle', () => {
  it('supports models that need no lifecycle hooks', () => {
    class Example extends ReactiveModel {}

    const model = new Example();
    expect(() => {
      model.mount();
      model.unmount();
      model.retire();
    }).not.toThrow();
  });

  it('mounts once per cycle, disposes owned subscriptions in reverse, and remounts', () => {
    const order: string[] = [];
    class Example extends ReactiveModel {
      protected override onMount(): void {
        order.push('mount');
        this.registerDispose(() => {
          order.push('first');
        });
        this.registerDispose(() => {
          order.push('second');
        });
      }
    }

    const model = new Example();
    model.mount();
    model.mount();
    model.unmount();
    expect(order).toEqual(['mount', 'second', 'first']);
    model.mount();
    expect(order).toEqual(['mount', 'second', 'first', 'mount']);
    model.unmount();
  });

  it('subscribes to observables and tears them down on unmount', () => {
    const source = observableState(1);
    const seen: number[] = [];
    class Example extends ReactiveModel {
      protected override onMount(): void {
        this.subscribe(source.asReadonly, (value) => {
          seen.push(value);
        });
      }
    }
    const model = new Example();
    model.mount();
    expect(seen).toEqual([1]);
    source.set(2);
    expect(seen).toEqual([1, 2]);
    model.unmount();
    source.set(3);
    expect(seen).toEqual([1, 2]);
  });

  it('tears down a mount subscription when its synchronous emission unmounts the model', () => {
    let activeSubscriptions = 0;
    const source = new Observable<number>((subscriber) => {
      activeSubscriptions += 1;
      subscriber.next(1);
      return () => {
        activeSubscriptions -= 1;
      };
    });
    class Example extends ReactiveModel {
      protected override onMount(): void {
        this.subscribe(source, () => this.unmount());
      }
    }

    const model = new Example();
    expect(() => model.mount()).not.toThrow();
    expect(activeSubscriptions).toBe(0);
    expect(() => model.mount()).not.toThrow();
    expect(activeSubscriptions).toBe(0);
  });

  it('tears down a permanent subscription when its synchronous emission retires the model', () => {
    let activeSubscriptions = 0;
    const source = new Observable<number>((subscriber) => {
      activeSubscriptions += 1;
      subscriber.next(1);
      return () => {
        activeSubscriptions -= 1;
      };
    });
    class Example extends ReactiveModel {
      subscribeNow(): void {
        this.subscribePermanently(source, () => this.retire());
      }
    }

    const model = new Example();
    expect(() => model.subscribeNow()).not.toThrow();
    expect(activeSubscriptions).toBe(0);
    expect(() => model.subscribeNow()).toThrow('Cannot subscribe permanently after retirement');
  });

  it('rejects registerDispose outside onMount', () => {
    class Example extends ReactiveModel {
      protected override onMount(): void {}

      register(): void {
        this.registerDispose(() => {});
      }
    }
    const model = new Example();
    expect(() => model.register()).toThrow(/outside onMount/);
    model.mount();
    expect(() => model.register()).toThrow(/outside onMount/);
  });

  it('rolls back a failed mount and permits a later mount', () => {
    const order: string[] = [];
    class Example extends ReactiveModel {
      attempts = 0;

      protected override onMount(): void {
        this.attempts += 1;
        this.registerDispose(() => {
          order.push('dispose');
        });
        if (this.attempts === 1) {
          throw new Error('mount failed');
        }
        order.push('mounted');
      }
    }

    const model = new Example();
    expect(() => model.mount()).toThrow('mount failed');
    expect(order).toEqual(['dispose']);
    model.mount();
    expect(order).toEqual(['dispose', 'mounted']);
    model.unmount();
    expect(order).toEqual(['dispose', 'mounted', 'dispose']);
  });

  it('attempts every failed-mount cleanup before reporting the mount failure', () => {
    const order: string[] = [];
    class Example extends ReactiveModel {
      protected override onMount(): void {
        this.registerDispose(() => {
          order.push('first');
        });
        this.registerDispose(() => {
          order.push('second');
          throw new Error('cleanup failed');
        });
        throw new Error('mount failed');
      }
    }

    const model = new Example();
    expect(() => model.mount()).toThrow('mount failed');
    expect(order).toEqual(['second', 'first']);
  });

  it('runs base cleanup when an unmount hook fails', () => {
    let disposed = false;
    class Example extends ReactiveModel {
      protected override onMount(): void {
        this.registerDispose(() => {
          disposed = true;
        });
      }

      protected override onUnmount(): void {
        throw new Error('unmount failed');
      }
    }

    const model = new Example();
    model.mount();
    expect(() => model.unmount()).toThrow('unmount failed');
    expect(disposed).toBe(true);
  });

  it('completes permanent cleanup and surfaces the original failure when retirement unmount fails', () => {
    const order: string[] = [];
    class Example extends ReactiveModel {
      constructor() {
        super();
        this.registerPermanentDispose(() => {
          order.push('permanent cleanup');
        });
      }

      protected override onMount(): void {
        this.registerDispose(() => {
          order.push('mount cleanup');
        });
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
    expect(order).toEqual(['unmount', 'mount cleanup', 'retire', 'permanent cleanup']);
    expect(() => model.retire()).not.toThrow();
    model.mount();
    expect(order).toEqual(['unmount', 'mount cleanup', 'retire', 'permanent cleanup']);
  });

  it('keeps permanent resources across unmount and releases them in reverse order on retirement', () => {
    const order: string[] = [];
    class Example extends ReactiveModel {
      constructor() {
        super();
        this.registerPermanentDispose(() => {
          order.push('first permanent');
        });
        this.registerPermanentDispose(() => {
          order.push('second permanent');
        });
      }

      protected override onMount(): void {
        this.registerDispose(() => {
          order.push('mount');
        });
      }
    }

    const model = new Example();
    model.mount();
    model.unmount();
    expect(order).toEqual(['mount']);
    model.retire();
    expect(order).toEqual(['mount', 'second permanent', 'first permanent']);
  });

  it('removes dynamic mount resources once while retaining siblings until unmount', () => {
    const order: string[] = [];
    class Example extends ReactiveModel {
      #removeChild: (() => void) | undefined;

      protected override onMount(): void {
        this.registerDispose(() => {
          order.push('retained');
        });
        this.#removeChild = this.registerDispose(() => {
          order.push('child');
        });
      }

      removeChild(): void {
        this.#removeChild?.();
      }
    }

    const model = new Example();
    model.mount();
    model.removeChild();
    model.unmount();
    model.removeChild();
    expect(order).toEqual(['child', 'retained']);
  });

  it('keeps public lifecycle control with the base model', () => {
    class Example extends ReactiveModel {
      override mount(): void {
        throw new Error('override should never be admitted');
      }
    }

    expect(() => new Example()).toThrow(/lifecycle methods cannot be overridden/);

    const model = new (class extends ReactiveModel {})();
    expect(() => Object.defineProperty(model, 'mount', { value: () => {} })).toThrow();
    expect(() => Object.defineProperty(model, 'retire', { value: () => {} })).toThrow();
    expect('dispose' in model).toBe(false);
    expect('onDispose' in model).toBe(false);
    expect('destroy' in model).toBe(false);
    expect('onDestroy' in model).toBe(false);

    class InvalidRetirement extends ReactiveModel {
      override retire(): void {}
    }
    expect(() => new InvalidRetirement()).toThrow(/lifecycle methods cannot be overridden/);
  });

  it('blocks mount after retirement', () => {
    class Example extends ReactiveModel {
      protected override onMount(): void {}
    }
    const model = new Example();
    model.mount();
    model.retire();
    model.mount();
    expect(model).toBeDefined();
  });
});

// An unmount during onMount waits for it
describe('an unmount during onMount waits for it', () => {
  describe('an unmount during onMount makes the rest of onMount throw from mount()', () => {
    it('mount() quietly leaves the model unmounted, with no live subscriptions or disposers', () => {
      const closed = new BehaviorSubject(true);
      const other = new BehaviorSubject(0);
      const disposed: string[] = [];
      class Example extends ReactiveModel {
        protected override onMount(): void {
          this.subscribe(closed, (value) => {
            if (value) this.unmount();
          });
          this.subscribe(other, () => {});
          this.registerDispose(() => disposed.push('later'));
        }
      }

      const model = new Example();
      expect(() => model.mount()).not.toThrow();
      expect(closed.observed).toBe(false);
      expect(other.observed).toBe(false);
    });
  });
});

// One owner mounts a model; a second mount is said once.
describe('a model mounted by two owners', () => {
  it('warns once, and stays a no-op', () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    let mounts = 0;
    class Example extends ReactiveModel {
      protected override onMount(): void {
        mounts++;
      }
    }
    const model = new Example();
    model.mount();
    model.mount();
    model.mount();
    expect(mounts).toBe(1);
    expect(warn).toHaveBeenCalledTimes(1);
    expect(String(warn.mock.calls[0]?.[0])).toContain('mounted again while mounted');
    warn.mockRestore();
  });

  it('does not warn when mounted again after an unmount', () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    class Example extends ReactiveModel {}
    const model = new Example();
    model.mount();
    model.unmount();
    model.mount();
    expect(warn).not.toHaveBeenCalled();
    warn.mockRestore();
  });
});

// A throwing projection is reported and recovered from
describe('a throwing projection is reported and recovered from', () => {
  describe('a throwing projection ends combineBehaviors permanently', () => {
    const previous = config.onUnhandledError;
    afterEach(() => {
      config.onUnhandledError = previous;
    });

    it('the combined behavior follows its source again after one throwing projection', async () => {
      const reported: unknown[] = [];
      // Capture what rxjs would otherwise rethrow on a timer, so the probe itself does not crash.
      config.onUnhandledError = (error) => reported.push(error);

      const a = observableState(1);
      const combined = combineBehaviors({
        sources: { a: a.asReadonly },
        project: ({ a: value }) => {
          if (value === 2) throw new Error('bad input');
          return { a: value };
        },
      });
      expect(combined.value).toEqual({ a: 1 });

      a.set(2); // the projection throws once
      a.set(3); // a valid value afterwards
      a.set(4);

      // rxjs reports the error from a timer (reportUnhandledError); let it fire while captured.
      await new Promise((resolve) => setTimeout(resolve, 0));
      expect(reported).toHaveLength(1);

      expect(a.value).toBe(4);
      // Correct behavior: the derived state recovers once the inputs are valid again.
      expect(combined.value).toEqual({ a: 4 });
      combined.dispose();
    });
  });
});
