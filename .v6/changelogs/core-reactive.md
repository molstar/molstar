## Core: reactive streams and states

- `@molstar/core/reactive`: a small rxjs-shaped stream library (`Observable`, `Subject`, `BehaviorSubject`,
  `Subscription`, `of`, `merge`, `combineLatest`, `timer`, `map`, `filter`, `distinctUntilChanged`, `skip`, `take`,
  `startWith`, `tap`, `switchMap`, `share`, `debounceTime`, `auditTime`, `throttleTime`, `firstValueFrom`,
  `lastValueFrom`) and states (`State`, `derive`, `batch`, `toState`) on one listener base. It interoperates with rxjs
  both ways (rxjs's `from`, `combineLatest`, `switchMap` take ours; ours take rxjs streams). One object per subscription;
  `derive` is glitch-free (a diamond computes and tells once) and lazy. Nothing uses it yet.
- `@molstar/core/reactive` also has `RecordState`, `derive` over named sources, `KeyedState` and `KeyedSelector` (a
  change visits only the keys it concerns), `Scope` and `ReactiveModel` (mount, unmount, retire; owned subscriptions,
  listeners and resources at one object each).
- `@molstar/core/util/object`: `isPlainRecord`, `sameShallow`, `sameData`, `shallowEqualObjects`,
  `shallowEqualArrays`, `deepEqual`; `shallowEqual` compares values by `Object.is` and makes no key arrays. The equality
  helpers of `@molstar/core/util` (`shallowEqual`, `shallowEqualObjects`, `shallowEqualArrays`, `deepEqual`)
  are the same implementations, with `Object.is` (`NaN` equals `NaN`, `-0` differs from `0`); its `arrayEqual` is
  `util/array`'s (`undefined` equals an empty array).
