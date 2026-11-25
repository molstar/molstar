/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as rx from 'rxjs';
import { describe, expect, expectTypeOf, it } from 'vitest';
import {
  BehaviorSubject,
  combineLatest,
  derive,
  from,
  isObservable,
  map,
  Observable,
  State,
  Subject,
  type Subscribable,
  type Subscription,
  switchMap,
  type Unsubscribable,
} from '../index.js';

// Interop with rxjs at the boundary (reactive §3.1): each direction, by type and at
// run time. Operators do not cross between the two `pipe`s; a mixed chain converts once with `from`.

describe('ours into rxjs', () => {
  it("rxjs's from, combineLatest and switchMap take ours (Symbol.observable)", () => {
    const s = new Subject<number>();
    const st = new State('a');
    const seen: unknown[] = [];
    const sub = rx
      .combineLatest([rx.from(s), st])
      .pipe(rx.switchMap(([n, t]) => new BehaviorSubject(`${t}${n}`)))
      .subscribe((v) => seen.push(v));
    s.next(1);
    st.set('b');
    expect(seen).toEqual(['a1', 'b1']);
    sub.unsubscribe();
    expect(s.observed).toBe(false);
    expect(st.observed).toBe(false);
  });

  it('ours is an rxjs Subscribable and InteropObservable; our Subject is an rxjs Observer', () => {
    expectTypeOf<Observable<number>>().toMatchTypeOf<rx.Subscribable<number>>();
    expectTypeOf<Subject<number>>().toMatchTypeOf<rx.Observer<number>>();
    expectTypeOf<Subscription>().toMatchTypeOf<rx.SubscriptionLike>();
    const target = new Subject<number>();
    const seen: number[] = [];
    target.subscribe((v) => seen.push(v));
    rx.of(1, 2).subscribe(target);
    expect(seen).toEqual([1, 2]);
    expect(target.isStopped).toBe(true);
  });

  it("rxjs's firstValueFrom takes ours through from", async () => {
    expect(await rx.firstValueFrom(rx.from(new BehaviorSubject(5)))).toBe(5);
  });
});

describe('rxjs into ours', () => {
  it('from, combineLatest, switchMap and derive take rxjs streams', () => {
    const rs = new rx.Subject<number>();
    const rb = new rx.BehaviorSubject('x');
    expectTypeOf<rx.Observable<number>>().toMatchTypeOf<Subscribable<number>>();
    const seen: unknown[] = [];
    const sub = combineLatest([from(rs), rb])
      .pipe(
        switchMap(([n, t]) => rx.of(`${t}${n}`)),
        map((v) => v.toUpperCase()),
      )
      .subscribe((v) => seen.push(v));
    rs.next(1);
    rb.next('y');
    expect(seen).toEqual(['X1', 'Y1']);
    sub.unsubscribe();
    expect(rs.observed).toBe(false);
  });

  it('isObservable is true for an rxjs stream too', () => {
    expect(isObservable(rx.of(1))).toBe(true);
    expect(isObservable(new Subject())).toBe(true);
    expect(isObservable({ subscribe() {} })).toBe(false);
  });

  it('derive follows a foreign state with a value and a subscribe (an rxjs BehaviorSubject)', () => {
    const rb = new rx.BehaviorSubject(1);
    const d = derive(rb as any as BehaviorSubject<number>, (x) => x * 2);
    const seen: number[] = [];
    const stop = d.listen((v) => seen.push(v));
    rb.next(2);
    stop();
    rb.next(3);
    expect(seen).toEqual([4]);
    expect(rb.observed).toBe(false);
    expect(d.value).toBe(6);
  });

  it('our subscribe returns rxjs-shaped subscriptions', () => {
    const u: Unsubscribable = new Subject<number>().subscribe(() => {});
    expectTypeOf(u).toMatchTypeOf<rx.Unsubscribable>();
  });
});
