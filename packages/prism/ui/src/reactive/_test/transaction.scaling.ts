/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { it } from 'vitest';
import { checkLinear, scaling } from '../../_test/scaling.js';
import { createOptimisticProjection, createUiTransactionCoordinator } from '../transaction.js';

// An optimistic projection with K contributions in flight grows linearly (`test:scaling`, ui/src/_test/scaling.ts): a
// contribution folds onto the projection as it stands, not the whole list again (a refold each: K² applies).

it('K contributions in flight on one projection, then retired at once', () => {
  const s = scaling(
    {
      setup: (n) => {
        const coordinator = createUiTransactionCoordinator();
        const projection = createOptimisticProjection(0, coordinator);
        const transactions = Array.from({ length: n }, () => coordinator.begin());
        return { coordinator, projection, transactions };
      },
      run: ({ coordinator, projection, transactions }) => {
        for (const t of transactions) projection.contribute(t, { apply: (v) => v + 1 });
        if (projection.projected.value !== transactions.length) throw new Error('fold');
        coordinator.reset();
      },
      teardown: ({ coordinator, projection }) => {
        projection.dispose();
        coordinator.dispose();
      },
    },
    2000,
  );
  checkLinear('optimistic projection: K contributions in flight', s);
});

it('K transactions retired while a projection holds K others', () => {
  // Every projection on a coordinator hears every retirement: one that holds nothing of it does nothing (no scan of its
  // contributions per event).
  const s = scaling(
    {
      setup: (n) => {
        const coordinator = createUiTransactionCoordinator();
        const projection = createOptimisticProjection(0, coordinator);
        for (let i = 0; i < n; i++) projection.contribute(coordinator.begin(), { apply: (v) => v + 1 });
        const elsewhere = Array.from({ length: n }, () => coordinator.begin());
        return { coordinator, projection, elsewhere };
      },
      run: ({ elsewhere }) => {
        for (const t of elsewhere) t.reject();
      },
      teardown: ({ coordinator, projection }) => {
        projection.dispose();
        coordinator.dispose();
      },
    },
    2000,
  );
  checkLinear('optimistic projection: K retired elsewhere, K held', s);
});
