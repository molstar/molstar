/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as rx from 'rxjs';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { ReactiveConfig } from '../config.js';
import * as ours from '../index.js';
import { BehaviorSubject, Observable, Subject, Subscription } from '../index.js';

// Each scenario is typed against rxjs; ours runs the same code (its API is the subset the scenarios use).
type Lib = typeof rx;
type Event = ['N', unknown] | ['E', unknown] | ['C'];

function record(source: { subscribe(o: any): { unsubscribe(): void } }) {
  const events: Event[] = [];
  const sub = source.subscribe({
    next: (v: unknown) => events.push(['N', v]),
    error: (e: unknown) => events.push(['E', e]),
    complete: () => events.push(['C']),
  });
  return { events, sub };
}

/** The same scenario through both libraries: rxjs is the oracle (reactive §3.2). */
function both(scenario: (lib: Lib) => unknown) {
  const expected = scenario(rx);
  const actual = scenario(ours as unknown as Lib);
  expect(actual).toEqual(expected);
  return actual;
}

afterEach(() => {
  ReactiveConfig.onUnhandledError = null;
  vi.useRealTimers();
});

describe('Subscription', () => {
  it('runs its teardowns once, in order, and runs one added after it ended at once', () => {
    const log: string[] = [];
    const s = new Subscription(() => log.push('initial'));
    s.add(() => log.push('a'));
    s.add({ unsubscribe: () => log.push('b') });
    s.unsubscribe();
    s.unsubscribe();
    s.add(() => log.push('late'));
    expect(log).toEqual(['initial', 'a', 'b', 'late']);
    expect(s.closed).toBe(true);
  });

  it('a child that ends first leaves its parent', () => {
    const parent = new Subscription();
    const log: string[] = [];
    for (let i = 0; i < 100; i++) {
      const child = new Subscription(() => log.push('child'));
      parent.add(child);
      child.unsubscribe();
    }
    expect((parent as any)._finalizers).toBeUndefined();
    parent.unsubscribe();
    expect(log.length).toBe(100);
  });

  it('runs every teardown when one throws, then throws', () => {
    const s = new Subscription();
    const log: string[] = [];
    s.add(() => {
      throw new Error('x');
    });
    s.add(() => log.push('after'));
    expect(() => s.unsubscribe()).toThrow('x');
    expect(log).toEqual(['after']);
  });

  it('EMPTY is closed', () => {
    expect(Subscription.EMPTY.closed).toBe(true);
  });
});

describe('Observable', () => {
  it('is cold: each subscribe runs the producer, and unsubscribe its teardown', () => {
    let runs = 0;
    let teardowns = 0;
    const o = new Observable<number>((s) => {
      runs++;
      s.next(runs);
      return () => teardowns++;
    });
    const a = record(o);
    const b = record(o);
    expect(a.events).toEqual([['N', 1]]);
    expect(b.events).toEqual([['N', 2]]);
    a.sub.unsubscribe();
    b.sub.unsubscribe();
    expect(teardowns).toBe(2);
  });

  it('ends at the first error or completion and tears down', () => {
    both((lib) => {
      const log: string[] = [];
      const o = new lib.Observable<number>((s) => {
        s.next(1);
        s.complete();
        s.next(2);
        s.error('late');
        return () => log.push('teardown');
      });
      return [record(o).events, log];
    });
  });

  it('a producer that throws is the stream error', () => {
    both(
      (lib) =>
        record(
          new lib.Observable<number>(() => {
            throw new Error('boom');
          }),
        ).events,
    );
  });

  it('reports what a subscriber throws, and goes on', () => {
    const reported: unknown[] = [];
    ReactiveConfig.onUnhandledError = (e) => reported.push(e);
    const s = new Subject<number>();
    const seen: number[] = [];
    s.subscribe(() => {
      throw new Error('bad');
    });
    s.subscribe((v) => seen.push(v));
    s.next(1);
    s.next(2);
    expect(reported).toEqual([new Error('bad'), new Error('bad')]);
    expect(seen).toEqual([1, 2]);
  });

  it('reports an error with no handler; rethrows it async when nothing is configured', () => {
    vi.useFakeTimers();
    const s = new Subject<number>();
    s.subscribe(() => {});
    s.error('lost');
    expect(() => vi.runAllTimers()).toThrow('lost');
  });

  it('pipes operators in order', () => {
    both((lib) => {
      const { map, filter } = lib;
      return record(
        lib.of(1, 2, 3, 4).pipe(
          map((x: number) => x * 10),
          filter((x: number) => x > 15),
          map((x: number, i: number) => `${i}:${x}`),
        ),
      ).events;
    });
  });
});

describe('Subject', () => {
  it('delivers nested: a value emitted by a listener reaches the others before the outer one', () => {
    both((lib) => {
      const s = new lib.Subject<number>();
      const seen: string[] = [];
      s.subscribe((v) => {
        seen.push(`a${v}`);
        if (v === 1) s.next(2);
      });
      s.subscribe((v) => seen.push(`b${v}`));
      s.next(1);
      return seen;
    });
  });

  it('a listener added during a delivery is not told that value; one removed before its turn is not told', () => {
    both((lib) => {
      const s = new lib.Subject<number>();
      const seen: string[] = [];
      let b: rx.Unsubscribable | undefined;
      s.subscribe((v) => {
        seen.push(`a${v}`);
        b?.unsubscribe();
        if (v === 1) s.subscribe((w) => seen.push(`c${w}`));
      });
      b = s.subscribe((v) => seen.push(`b${v}`));
      s.next(1);
      s.next(2);
      return seen;
    });
  });

  it('ends: later next ignored, a late subscriber gets the end at once', () => {
    both((lib) => {
      const s = new lib.Subject<number>();
      const early = record(s);
      s.next(1);
      s.complete();
      s.next(2);
      const late = record(s);
      const e = new lib.Subject<number>();
      e.error('x');
      return [early.events, late.events, record(e).events, s.observed];
    });
  });

  it('observed, and asObservable hides next', () => {
    const s = new Subject<number>();
    expect(s.observed).toBe(false);
    const view = s.asObservable();
    expect('next' in view).toBe(false);
    const r = record(view);
    expect(s.observed).toBe(true);
    s.next(1);
    r.sub.unsubscribe();
    expect(s.observed).toBe(false);
    expect(r.events).toEqual([['N', 1]]);
  });

  it('listen is a plain callback until stopped', () => {
    const s = new Subject<number>();
    const seen: number[] = [];
    const stop = s.listen((v) => seen.push(v));
    s.next(1);
    stop();
    s.next(2);
    expect(seen).toEqual([1]);
    expect(s.observed).toBe(false);
  });

  it('keeps its list compact when many listeners go (holes dropped at half)', () => {
    const s = new Subject<number>();
    const subs = Array.from({ length: 1000 }, () => s.subscribe(() => {}));
    for (let i = 0; i < 900; i++) subs[i].unsubscribe();
    expect((s as any)._subs.length).toBeLessThan(400);
    let told = 0;
    for (let i = 900; i < 1000; i++) subs[i].unsubscribe();
    s.subscribe(() => told++);
    s.next(1);
    expect(told).toBe(1);
    expect((s as any)._subs.length).toBe(1);
  });

  it('unsubscribe closes it and drops its listeners', () => {
    const s = new Subject<number>();
    const r = record(s);
    s.unsubscribe();
    s.next(1);
    expect(r.events).toEqual([]);
    expect(s.closed).toBe(true);
    expect(s.observed).toBe(false);
  });
});

describe('BehaviorSubject', () => {
  it('replays its value, nested, then follows', () => {
    both((lib) => {
      const s = new lib.BehaviorSubject(1);
      const seen: string[] = [];
      s.subscribe((v) => {
        seen.push(`a${v}`);
        if (v === 1) s.next(2);
      });
      s.subscribe((v) => seen.push(`b${v}`));
      s.next(3);
      return [seen, s.value, s.getValue()];
    });
  });

  it('getValue throws after an error; a completed one replays nothing', () => {
    both((lib) => {
      const s = new lib.BehaviorSubject(1);
      s.complete();
      const r = record(s);
      const e = new lib.BehaviorSubject(1);
      e.error('x');
      let thrown: unknown;
      try {
        e.getValue();
      } catch (err) {
        thrown = err;
      }
      return [r.events, s.getValue(), thrown];
    });
  });

  it('is a ReadonlyState: listen does not replay', () => {
    const s = new BehaviorSubject(1);
    const seen: number[] = [];
    s.listen((v) => seen.push(v));
    s.next(2);
    expect(seen).toEqual([2]);
    expect(ours.isState(s)).toBe(true);
  });
});

describe('creation', () => {
  it('of, EMPTY, from (array, iterable, observable)', () => {
    both((lib) => [
      record(lib.of(1, 2)).events,
      record(lib.EMPTY).events,
      record(lib.from([3, 4])).events,
      record(lib.from(new Set([5]))).events,
      record(lib.from(lib.of(6))).events,
    ]);
  });

  it('from a promise', async () => {
    const r = record(ours.from(Promise.resolve(7)));
    await Promise.resolve();
    await Promise.resolve();
    expect(r.events).toEqual([['N', 7], ['C']]);
  });

  it('merge: every value as it comes, completes when all have, errors with the first', () => {
    both((lib) => {
      const a = new lib.Subject<number>();
      const b = new lib.Subject<number>();
      const r = record(lib.merge(a, b));
      a.next(1);
      b.next(2);
      a.complete();
      b.next(3);
      const before = r.events.slice();
      b.complete();
      const c = new lib.Subject<number>();
      const e = record(lib.merge(c, lib.of(9)));
      c.error('x');
      return [before, r.events, e.events, record(lib.merge()).events];
    });
  });

  it('merge unsubscribes every source', () => {
    const a = new Subject<number>();
    const b = new Subject<number>();
    const sub = ours.merge(a, b).subscribe(() => {});
    expect(a.observed && b.observed).toBe(true);
    sub.unsubscribe();
    expect(a.observed || b.observed).toBe(false);
  });

  it('combineLatest: from when all told, a new array each time, completes when all have', () => {
    both((lib) => {
      const a = new lib.Subject<number>();
      const b = new lib.BehaviorSubject('x');
      const arrays: unknown[] = [];
      const events: Event[] = [];
      lib.combineLatest([a, b]).subscribe({
        next: (v) => {
          arrays.push(v);
          events.push(['N', v.slice()]);
        },
        complete: () => events.push(['C']),
      });
      a.next(1);
      b.next('y');
      a.next(2);
      a.complete();
      b.complete();
      return [events, new Set(arrays).size === arrays.length, record(lib.combineLatest([])).events];
    });
  });

  it('timer: once, or with a period', () => {
    vi.useFakeTimers();
    both((lib) => {
      const once = record(lib.timer(10));
      vi.advanceTimersByTime(9);
      const early = once.events.slice();
      vi.advanceTimersByTime(1);
      const periodic = record(lib.timer(5, 3));
      vi.advanceTimersByTime(12);
      periodic.sub.unsubscribe();
      vi.advanceTimersByTime(100);
      return [early, once.events, periodic.events];
    });
  });
});

describe('operators', () => {
  it('map and filter pass the index; a throwing callback is the stream error', () => {
    both((lib) => [
      record(lib.of(5, 6, 7).pipe(lib.map((x: number, i: number) => x + i))).events,
      record(lib.of(5, 6, 7).pipe(lib.filter((_: number, i: number) => i !== 1))).events,
      record(
        lib.of(1, 2).pipe(
          lib.map((x: number) => {
            if (x === 2) throw new Error('bad');
            return x;
          }),
        ),
      ).events,
      record(lib.of(0, 1, '', 'a', null).pipe(lib.filter(Boolean))).events,
    ]);
  });

  it('distinctUntilChanged: ===, a comparator, a key selector', () => {
    both((lib) => [
      record(lib.of(1, 1, 2, 2, 1).pipe(lib.distinctUntilChanged())).events,
      record(lib.of(1, 2, 4, 5).pipe(lib.distinctUntilChanged((a: number, b: number) => Math.abs(a - b) < 2))).events,
      record(
        lib.of({ k: 1 }, { k: 1 }, { k: 2 }).pipe(
          lib.distinctUntilChanged(
            (a: number, b: number) => a === b,
            (o: { k: number }) => o.k,
          ),
        ),
      ).events,
    ]);
  });

  it('skip, take, startWith', () => {
    both((lib) => {
      const s = new lib.Subject<number>();
      const taken = record(s.pipe(lib.take(2)));
      s.next(1);
      s.next(2);
      const observedAfterTake = s.observed;
      return [
        record(lib.of(1, 2, 3).pipe(lib.skip(1))).events,
        record(lib.of(1, 2, 3).pipe(lib.take(0))).events,
        taken.events,
        observedAfterTake,
        record(lib.of(3).pipe(lib.startWith(1, 2))).events,
      ];
    });
  });

  it('tap sees each notification; a throwing tap is the stream error', () => {
    both((lib) => {
      const log: unknown[] = [];
      const ok = record(
        lib.of(1, 2).pipe(lib.tap({ next: (v: number) => log.push(v), complete: () => log.push('done') })),
      ).events;
      const bad = record(
        lib.of(1).pipe(
          lib.tap(() => {
            throw new Error('x');
          }),
        ),
      ).events;
      return [ok, log, bad];
    });
  });

  it('switchMap follows only the latest inner; completes when the source and the inner have', () => {
    both((lib) => {
      const outer = new lib.Subject<number>();
      const inners = [new lib.Subject<string>(), new lib.Subject<string>()];
      const r = record(outer.pipe(lib.switchMap((i: number) => inners[i])));
      outer.next(0);
      inners[0].next('a');
      outer.next(1);
      inners[0].next('stale');
      inners[1].next('b');
      outer.complete();
      const beforeInnerDone = r.events.slice();
      inners[1].complete();
      return [beforeInnerDone, r.events, inners[0].observed, inners[1].observed];
    });
  });

  it('switchMap takes promises, arrays and of', async () => {
    both((lib) => record(lib.of(1, 2).pipe(lib.switchMap((x: number) => [x, x * 10]))).events);
    const r = record(ours.of(1).pipe(ours.switchMap((x) => Promise.resolve(x + 1))));
    await new Promise((resolve) => setTimeout(resolve, 0));
    expect(r.events).toEqual([['N', 2], ['C']]);
  });

  it('share: one upstream subscription, made by the first and ended by the last, reset after the end', () => {
    let subscribes = 0;
    let teardowns = 0;
    const source = new Subject<number>();
    const counted = new Observable<number>((s) => {
      subscribes++;
      const sub = source.subscribe(s);
      return () => {
        teardowns++;
        sub.unsubscribe();
      };
    });
    const shared = counted.pipe(ours.share());
    const a = record(shared);
    const b = record(shared);
    source.next(1);
    expect(subscribes).toBe(1);
    a.sub.unsubscribe();
    b.sub.unsubscribe();
    expect(teardowns).toBe(1);
    record(shared);
    expect(subscribes).toBe(2);
    expect(a.events).toEqual([['N', 1]]);
  });

  it('share passes a synchronous source and its end to the first subscriber, as rxjs', () => {
    both((lib) => {
      const shared = lib.of(1, 2).pipe(lib.share());
      return [record(shared).events, record(shared).events];
    });
  });
});

describe('time operators (rxjs as the oracle, fake timers)', () => {
  type Step = [number, number | 'C'];

  /** Feeds `steps` (at their times) through `op` in each library and records what comes out, with the time. */
  function timed(steps: Step[], op: (lib: Lib) => rx.MonoTypeOperatorFunction<any>, until = 1000) {
    return both((lib) => {
      vi.useFakeTimers();
      const start = Date.now();
      const s = new lib.Subject<number>();
      const out: unknown[] = [];
      (s.pipe(op(lib) as any) as rx.Observable<number>).subscribe({
        next: (v) => out.push([Date.now() - start, v]),
        complete: () => out.push([Date.now() - start, 'C']),
      });
      let t = 0;
      for (const [at, v] of steps) {
        vi.advanceTimersByTime(at - t);
        t = at;
        if (v === 'C') s.complete();
        else s.next(v);
      }
      vi.advanceTimersByTime(until - t);
      vi.useRealTimers();
      return out;
    });
  }

  const burst: Step[] = [
    [0, 1],
    [10, 2],
    [20, 3],
    [120, 4],
    [130, 5],
    [300, 6],
  ];

  it('debounceTime', () => {
    timed(burst, (lib) => lib.debounceTime(50));
    timed([...burst, [310, 'C']], (lib) => lib.debounceTime(50));
  });

  it('auditTime', () => {
    timed(burst, (lib) => lib.auditTime(50));
    timed([...burst, [310, 'C']], (lib) => lib.auditTime(50));
  });

  for (const config of [
    { leading: true, trailing: true },
    { leading: false, trailing: true },
    { leading: true, trailing: false },
    undefined,
  ]) {
    it(`throttleTime ${JSON.stringify(config)}`, () => {
      timed(burst, (lib) => lib.throttleTime(50, undefined, config));
      timed([...burst, [310, 'C']], (lib) => lib.throttleTime(50, undefined, config));
      timed([...burst.slice(0, 5), [140, 'C']], (lib) => lib.throttleTime(50, undefined, config));
    });
  }

  it('random sequences agree with rxjs', () => {
    let seed = 7;
    const rand = () => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed / 0x7fffffff;
    };
    for (let round = 0; round < 30; round++) {
      const steps: Step[] = [];
      let t = 0;
      for (let i = 0; i < 12; i++) {
        t += Math.floor(rand() * 40);
        steps.push([t, i]);
      }
      if (rand() < 0.5) steps.push([t + Math.floor(rand() * 40), 'C']);
      timed(steps, (lib) => lib.debounceTime(25));
      timed(steps, (lib) => lib.auditTime(25));
      timed(steps, (lib) => lib.throttleTime(25, undefined, { leading: true, trailing: true }));
      timed(steps, (lib) => lib.throttleTime(25, undefined, { leading: false, trailing: true }));
    }
  });

  it('a timer is cleared on unsubscribe', () => {
    vi.useFakeTimers();
    const s = new Subject<number>();
    const seen: number[] = [];
    const sub = s.pipe(ours.debounceTime(10)).subscribe((v) => seen.push(v));
    s.next(1);
    sub.unsubscribe();
    expect(vi.getTimerCount()).toBe(0);
    vi.advanceTimersByTime(20);
    expect(seen).toEqual([]);
  });
});

describe('promises', () => {
  it('firstValueFrom and lastValueFrom', async () => {
    const s = new Subject<number>();
    const first = ours.firstValueFrom(s);
    s.next(1);
    expect(await first).toBe(1);
    expect(s.observed).toBe(false);
    expect(await ours.lastValueFrom(ours.of(1, 2, 3))).toBe(3);
    await expect(ours.firstValueFrom(ours.EMPTY)).rejects.toBeInstanceOf(ours.EmptyError);
    expect(await ours.firstValueFrom(ours.EMPTY, { defaultValue: 'd' })).toBe('d');
    await expect(ours.lastValueFrom(ours.from(Promise.reject('x')))).rejects.toBe('x');
  });
});

describe('review fixes', () => {
  it('Subject.unsubscribe drops every listener (no skips while the list compacts)', () => {
    const s = new Subject<number>();
    for (let i = 0; i < 7; i++) s.subscribe(() => {});
    s.unsubscribe();
    expect(s.observed).toBe(false);
    expect((s as any)._subs.length).toBe(0);
  });

  it('a time operator that completes after its downstream ended does not throw', () => {
    const reported: unknown[] = [];
    ReactiveConfig.onUnhandledError = (e) => reported.push(e);
    const seen: number[] = [];
    ours
      .of(1)
      .pipe(ours.debounceTime(10), ours.take(1))
      .subscribe((v) => seen.push(v));
    expect(seen).toEqual([1]);
    expect(reported).toEqual([]);
  });

  it('asObservable is made once', () => {
    const s = new Subject<number>();
    expect(s.asObservable()).toBe(s.asObservable());
  });
});
