/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { createOptimisticProjection, createUiTransactionCoordinator, type UiTransactionEvent } from '../transaction.js';

describe('UI transaction coordinator', () => {
  it('begins isolated transactions and rejects new work after disposal', () => {
    const coordinator = createUiTransactionCoordinator();
    expect(coordinator.begin().id).toBe('tx-1');
    expect(coordinator.begin().id).toBe('tx-2');

    coordinator.dispose();

    expect(() => coordinator.begin()).toThrow('UI transaction coordinator is disposed');
  });

  it('reports acknowledgement without ending a transaction', () => {
    const coordinator = createUiTransactionCoordinator();
    const events: UiTransactionEvent[] = [];
    coordinator.events.subscribe((event) => events.push(event));
    const transaction = coordinator.begin();

    transaction.acknowledge();
    transaction.acknowledge();
    transaction.reject();

    expect(events).toEqual([
      { type: 'acknowledged', id: transaction.id },
      { type: 'rejected', id: transaction.id },
    ]);
  });
});

describe('optimistic projection', () => {
  it('folds overlapping contributions in stable order over the latest authoritative base', () => {
    const coordinator = createUiTransactionCoordinator();
    const projection = createOptimisticProjection(['base'], coordinator);
    const first = coordinator.begin();
    const second = coordinator.begin();

    projection.contribute(first, { apply: (value) => [...value, 'first'] });
    projection.contribute(second, { apply: (value) => [...value, 'second'] });
    projection.setBase(['newer']);

    expect(projection.projected.value).toEqual(['newer', 'first', 'second']);
  });

  it('preserves a contribution after receipt until authoritative absorption', () => {
    const coordinator = createUiTransactionCoordinator();
    const projection = createOptimisticProjection('base', coordinator);
    const events: UiTransactionEvent[] = [];
    coordinator.events.subscribe((event) => events.push(event));
    const transaction = coordinator.begin();
    projection.contribute(transaction, {
      apply: (value) => `${value}+one`,
      isAbsorbed: (value) => value.includes('one'),
    });

    transaction.acknowledge();
    expect(projection.projected.value).toBe('base+one');

    projection.setBase('base+one');
    expect(projection.projected.value).toBe('base+one');
    expect(events).toEqual([
      { type: 'acknowledged', id: transaction.id },
      { type: 'absorbed', id: transaction.id },
    ]);
  });

  it('accepts authoritative absorption before a reordered acknowledgement', () => {
    const coordinator = createUiTransactionCoordinator();
    const projection = createOptimisticProjection('base', coordinator);
    const events: UiTransactionEvent[] = [];
    coordinator.events.subscribe((event) => events.push(event));
    const transaction = coordinator.begin();
    projection.contribute(transaction, {
      apply: (value) => `${value}+one`,
      isAbsorbed: (value) => value.endsWith('+one'),
    });

    projection.setBase('base+one');
    transaction.acknowledge();

    expect(projection.projected.value).toBe('base+one');
    expect(events).toEqual([{ type: 'absorbed', id: transaction.id }]);
  });

  it('ends a transaction only after all of its contributions are absorbed', () => {
    const coordinator = createUiTransactionCoordinator();
    const projection = createOptimisticProjection<string[]>([], coordinator);
    const events: UiTransactionEvent[] = [];
    coordinator.events.subscribe((event) => events.push(event));
    const transaction = coordinator.begin();
    projection.contribute(transaction, {
      apply: (value) => [...value, 'first'],
      isAbsorbed: (value) => value.includes('first'),
    });
    projection.contribute(transaction, {
      apply: (value) => [...value, 'second'],
      isAbsorbed: (value) => value.includes('second'),
    });

    projection.setBase(['first']);
    expect(projection.projected.value).toEqual(['first', 'second']);
    expect(events).toEqual([]);

    projection.setBase(['first', 'second']);
    expect(projection.projected.value).toEqual(['first', 'second']);
    expect(events).toEqual([{ type: 'absorbed', id: transaction.id }]);
  });

  it('reconciles after any attachment is released without absorption', () => {
    const coordinator = createUiTransactionCoordinator();
    const firstProjection = createOptimisticProjection('first', coordinator);
    const secondProjection = createOptimisticProjection('second', coordinator);
    const events: UiTransactionEvent[] = [];
    coordinator.events.subscribe((event) => events.push(event));
    const transaction = coordinator.begin();
    firstProjection.contribute(transaction, { apply: (value) => `${value}+pending` });
    secondProjection.contribute(transaction, {
      apply: (value) => `${value}+pending`,
      isAbsorbed: (value) => value.endsWith('+pending'),
    });

    firstProjection.reset();
    expect(events).toEqual([]);
    secondProjection.setBase('second+pending');

    expect(events).toEqual([{ type: 'reconciled', id: transaction.id }]);
  });

  it('rejects one transaction without dropping newer work or rewinding the base', () => {
    const coordinator = createUiTransactionCoordinator();
    const projection = createOptimisticProjection(['base'], coordinator);
    const first = coordinator.begin();
    const second = coordinator.begin();
    projection.contribute(first, { apply: (value) => [...value, 'first'] });
    projection.contribute(second, { apply: (value) => [...value, 'second'] });
    projection.setBase(['base', 'remote']);

    first.reject();

    expect(projection.projected.value).toEqual(['base', 'remote', 'second']);
  });

  it('reconciles a divergent authoritative result without disturbing another transaction', () => {
    const coordinator = createUiTransactionCoordinator();
    const projection = createOptimisticProjection(['base'], coordinator);
    const first = coordinator.begin();
    const second = coordinator.begin();
    projection.contribute(first, { apply: (value) => [...value, 'predicted'] });
    projection.contribute(second, { apply: (value) => [...value, 'later'] });

    projection.setBase(['base', 'authoritative']);
    first.reconcile();

    expect(projection.projected.value).toEqual(['base', 'authoritative', 'later']);
  });

  it('reports an inapplicable contribution and retains it for owner resolution', () => {
    const coordinator = createUiTransactionCoordinator();
    const projection = createOptimisticProjection({ version: 1 }, coordinator);
    const divergence = vi.fn();
    projection.divergences.subscribe(divergence);
    const transaction = coordinator.begin();
    projection.contribute(transaction, {
      apply: (base) => {
        if (base.version !== 1) throw new Error('unsupported version');
        return { version: 2 };
      },
    });

    projection.setBase({ version: 3 });

    expect(projection.projected.value).toEqual({ version: 3 });
    expect(divergence).toHaveBeenCalledWith({
      transactionId: transaction.id,
      phase: 'application',
      authoritativeBase: { version: 3 },
      cause: expect.objectContaining({ message: 'unsupported version' }),
    });

    transaction.reconcile();
    expect(projection.projected.value).toEqual({ version: 3 });
  });

  it('resets every projection to its latest authoritative base', () => {
    const coordinator = createUiTransactionCoordinator();
    const firstProjection = createOptimisticProjection('first', coordinator);
    const secondProjection = createOptimisticProjection('second', coordinator);
    const first = coordinator.begin();
    const second = coordinator.begin();
    firstProjection.contribute(first, { apply: (value) => `${value}+pending` });
    secondProjection.contribute(second, { apply: (value) => `${value}+pending` });
    firstProjection.setBase('first-newer');
    secondProjection.setBase('second-newer');

    coordinator.reset();

    expect(firstProjection.projected.value).toBe('first-newer');
    expect(secondProjection.projected.value).toBe('second-newer');
    expect(() => firstProjection.contribute(first, { apply: (value) => value })).toThrow(
      'UI transaction is no longer active',
    );
  });

  it('resets one projection without disturbing another projection on the same coordinator', () => {
    const coordinator = createUiTransactionCoordinator();
    const firstProjection = createOptimisticProjection('first', coordinator);
    const secondProjection = createOptimisticProjection('second', coordinator);
    const first = coordinator.begin();
    const second = coordinator.begin();
    firstProjection.contribute(first, { apply: (value) => `${value}+pending` });
    secondProjection.contribute(second, { apply: (value) => `${value}+pending` });

    firstProjection.reset();

    expect(firstProjection.projected.value).toBe('first');
    expect(secondProjection.projected.value).toBe('second+pending');
    expect(() => firstProjection.contribute(first, { apply: (value) => value })).toThrow(
      'UI transaction is no longer active',
    );
    second.reject();
  });

  it('enforces coordinator ownership and deterministic late attachment', () => {
    const owner = createUiTransactionCoordinator();
    const other = createUiTransactionCoordinator();
    const projection = createOptimisticProjection('base', owner);
    const foreign = other.begin();
    const ended = owner.begin();
    ended.reject();

    expect(() => projection.contribute(foreign, { apply: (value) => value })).toThrow(
      'UI transaction belongs to a different coordinator',
    );
    expect(() => projection.contribute(ended, { apply: (value) => value })).toThrow(
      'UI transaction is no longer active',
    );
  });

  it('disposes transient state once and ignores later base updates', () => {
    const coordinator = createUiTransactionCoordinator();
    const projection = createOptimisticProjection('base', coordinator);
    const transaction = coordinator.begin();
    projection.contribute(transaction, { apply: (value) => `${value}+pending` });
    projection.setBase('newer');

    coordinator.dispose();
    projection.setBase('late');
    coordinator.dispose();

    expect(projection.projected.value).toBe('newer');
    expect(() => projection.contribute(transaction, { apply: (value) => value })).toThrow(
      'Optimistic projection is disposed',
    );
  });
});

// An inapplicable contribution hides no later one
describe('an inapplicable contribution hides no later one', () => {
  function setup() {
    const coordinator = createUiTransactionCoordinator();
    const projection = createOptimisticProjection<number[]>([], coordinator);
    const first = coordinator.begin();
    const second = coordinator.begin();
    projection.contribute(first, {
      apply: (base) => {
        if (base.length > 0) throw new Error('no longer applies');
        return [...base, 1];
      },
    });
    projection.contribute(second, { apply: (base) => [...base, 2] });
    return { coordinator, projection, first, second };
  }

  describe('optimistic projection: an inapplicable contribution', () => {
    it('does not hide a later, applicable contribution while the divergence is unresolved', () => {
      const { projection } = setup();
      expect(projection.projected.value).toEqual([1, 2]);

      projection.setBase([7]);

      // The later transaction is still active and applies to [7]; it must stay visible.
      expect(projection.projected.value).toEqual([7, 2]);
    });

    it('shows the later contribution on every newer base while the divergence stays unresolved', () => {
      const { projection } = setup();
      projection.setBase([7]);
      projection.setBase([8]);
      projection.setBase([9]);

      expect(projection.projected.value).toEqual([9, 2]);
    });

    it('does not deliver a partial projection to subscribers even when the owner resolves synchronously', () => {
      const { projection, first } = setup();
      const seen: number[][] = [];
      projection.projected.subscribe((v) => seen.push(v));
      // The owner rejects the diverging transaction at once.
      projection.divergences.subscribe(() => first.reject());

      projection.setBase([7]);

      // Final state is right ...
      expect(projection.projected.value).toEqual([7, 2]);
      // ... but subscribers must never have seen [7] (the base without the still-active second edit).
      expect(seen).not.toContainEqual([7]);
    });
  });
});

// setBase keeps changes made by divergence handlers
describe('setBase keeps changes made by divergence handlers', () => {
  describe('setBase overwrites contributions changed by a re-entrant divergence subscriber', () => {
    it('a contribution rejected from an absorption divergence stays removed', () => {
      const coordinator = createUiTransactionCoordinator();
      const projection = createOptimisticProjection<readonly number[]>([0], coordinator);
      const events: UiTransactionEvent[] = [];
      coordinator.events.subscribe((e) => events.push(e));
      const tx = coordinator.begin();
      projection.contribute(tx, {
        apply: (base) => [...base, 99],
        isAbsorbed: () => {
          throw new Error('cannot tell');
        },
      });
      // The owner resolves the divergence by rejecting the edit.
      projection.divergences.subscribe((d) => {
        if (d.transactionId === tx.id) tx.reject();
      });

      projection.setBase([1]);

      expect(events).toContainEqual({ type: 'rejected', id: tx.id });
      // Correct: the rejected edit is no longer projected. Defect: [1, 99].
      expect(projection.projected.value).toEqual([1]);
      projection.setBase([2]);
      expect(projection.projected.value).toEqual([2]);
    });

    it('a contribution added from an absorption divergence is kept and retires its transaction', () => {
      const coordinator = createUiTransactionCoordinator();
      const projection = createOptimisticProjection<readonly number[]>([0], coordinator);
      const events: UiTransactionEvent[] = [];
      coordinator.events.subscribe((e) => events.push(e));
      const tx1 = coordinator.begin();
      projection.contribute(tx1, {
        apply: (base) => base,
        isAbsorbed: () => {
          throw new Error('cannot tell');
        },
      });
      const tx2 = coordinator.begin();
      let added = false;
      projection.divergences.subscribe(() => {
        if (added) return;
        added = true;
        projection.contribute(tx2, {
          apply: (base) => [...base, 7],
          isAbsorbed: (base) => base.includes(7),
        });
      });

      projection.setBase([1]);
      // Correct: tx2's edit is projected. Defect: [1].
      expect(projection.projected.value).toEqual([1, 7]);

      // Once the authoritative base absorbs it, tx2 retires. Defect: no event ever comes.
      projection.setBase([1, 7]);
      expect(events).toContainEqual({ type: 'absorbed', id: tx2.id });
    });
  });
});

describe('an appended contribution folds onto the projection as it stands', () => {
  it('gives what a refold gives, and reports a contribution still failing again', () => {
    const coordinator = createUiTransactionCoordinator();
    const projection = createOptimisticProjection<number[]>([], coordinator);
    const seen: number[][] = [];
    projection.projected.subscribe((v) => seen.push(v));
    const reported: string[] = [];
    projection.divergences.subscribe((d) => reported.push(d.transactionId));
    const applied: number[] = [];
    const failing = coordinator.begin();
    projection.contribute(failing, {
      apply: () => {
        throw new Error('never applies');
      },
    });
    const txs = [1, 2, 3].map((k) => {
      const t = coordinator.begin();
      projection.contribute(t, {
        apply: (v) => {
          applied.push(k);
          return [...v, k];
        },
      });
      return t;
    });
    // Each contribution applied once, at its arrival; the failing one reported at each publish.
    expect(applied).toEqual([1, 2, 3]);
    expect(projection.projected.value).toEqual([1, 2, 3]);
    expect(reported).toEqual([failing.id, failing.id, failing.id, failing.id]);
    // A retirement refolds what is left over the base.
    txs[1].reject();
    expect(projection.projected.value).toEqual([1, 3]);
    expect(seen.at(-1)).toEqual([1, 3]);
    failing.reject();
    expect(projection.projected.value).toEqual([1, 3]);
  });

  it('a retirement of a transaction it holds nothing of publishes nothing', () => {
    const coordinator = createUiTransactionCoordinator();
    const projection = createOptimisticProjection(0, coordinator);
    projection.contribute(coordinator.begin(), { apply: (v) => v + 1 });
    let published = 0;
    projection.projected.subscribe(() => published++);
    const before = published;
    coordinator.begin().reject();
    expect(published).toBe(before);
    expect(projection.projected.value).toBe(1);
  });
});
