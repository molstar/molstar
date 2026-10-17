/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Observable, Subject } from 'rxjs';

import type { Disposable } from './disposable-group.js';
import { observableState } from './observable-state.js';
import type { ReadonlyState } from './readonly-state.js';

declare const uiTransactionIdBrand: unique symbol;

export type UiTransactionId = string & { readonly [uiTransactionIdBrand]: undefined };

export type UiTransactionEvent =
  | {
      readonly type: 'acknowledged' | 'absorbed' | 'rejected' | 'reconciled';
      readonly id: UiTransactionId;
    }
  | { readonly type: 'reset' | 'disposed' };

export type UiTransaction = {
  readonly id: UiTransactionId;
  /** Records request receipt without resolving this transaction or removing its contributions. */
  acknowledge(): void;
  /** Retires this transaction after the owner rejects its predicted outcome. */
  reject(): void;
  /** Retires this transaction after the owner explicitly reconciles an authoritative result. */
  reconcile(): void;
};

export type UiTransactionCoordinator = Disposable & {
  readonly events: Observable<UiTransactionEvent>;
  begin(): UiTransaction;
  reset(): void;
};

export type OptimisticContribution<T> = {
  apply(base: T): T;
  isAbsorbed?(base: T): boolean;
};

export type OptimisticProjectionDivergence<T> = {
  readonly transactionId: UiTransactionId;
  readonly phase: 'absorption' | 'application';
  readonly authoritativeBase: T;
  readonly cause: unknown;
};

export type OptimisticProjection<T> = Disposable & {
  readonly projected: ReadonlyState<T>;
  readonly divergences: Observable<OptimisticProjectionDivergence<T>>;
  setBase(base: T): void;
  contribute(transaction: UiTransaction, contribution: OptimisticContribution<T>): void;
  reset(): void;
};

type TransactionState = {
  readonly coordinator: CoordinatorState;
  readonly id: UiTransactionId;
  active: boolean;
  acknowledged: boolean;
  attachments: number;
  releasedUnabsorbed: boolean;
};

type CoordinatorState = {
  readonly events: Subject<UiTransactionEvent>;
  readonly transactions: Set<TransactionState>;
  disposed: boolean;
  sequence: number;
};

type TrackedContribution<T> = OptimisticContribution<T> & {
  readonly transactionId: UiTransactionId;
  release(absorbed: boolean): void;
  /** Out of the projection (absorbed, dropped, cleared): reported no more. */
  dropped: boolean;
};

type Failure<T> = { readonly contribution: TrackedContribution<T>; readonly cause: unknown };

const transactionStates = new WeakMap<UiTransaction, TransactionState>();
const coordinatorStates = new WeakMap<UiTransactionCoordinator, CoordinatorState>();

function uiTransactionId(value: string): UiTransactionId {
  return value as UiTransactionId;
}

function finishTransaction(state: TransactionState, type: 'absorbed' | 'rejected' | 'reconciled'): void {
  if (!state.active) {
    return;
  }
  state.active = false;
  state.attachments = 0;
  state.coordinator.transactions.delete(state);
  state.coordinator.events.next({ type, id: state.id });
}

function attachTransaction(
  coordinator: UiTransactionCoordinator,
  transaction: UiTransaction,
): { readonly id: UiTransactionId; release(absorbed: boolean): void } {
  const coordinatorState = coordinatorStates.get(coordinator);
  const transactionState = transactionStates.get(transaction);
  if (coordinatorState === undefined || coordinatorState.disposed) {
    throw new Error('UI transaction coordinator is disposed');
  }
  if (transactionState === undefined || transactionState.coordinator !== coordinatorState) {
    throw new Error('UI transaction belongs to a different coordinator');
  }
  if (!transactionState.active) {
    throw new Error('UI transaction is no longer active');
  }

  transactionState.attachments += 1;
  let attached = true;
  return {
    id: transactionState.id,
    release(absorbed) {
      if (!attached) {
        return;
      }
      attached = false;
      if (!transactionState.active) {
        return;
      }
      if (!absorbed) {
        transactionState.releasedUnabsorbed = true;
      }
      transactionState.attachments -= 1;
      if (transactionState.attachments === 0) {
        finishTransaction(transactionState, transactionState.releasedUnabsorbed ? 'reconciled' : 'absorbed');
      }
    },
  };
}

export function createUiTransactionCoordinator(): UiTransactionCoordinator {
  const state: CoordinatorState = {
    events: new Subject<UiTransactionEvent>(),
    transactions: new Set<TransactionState>(),
    disposed: false,
    sequence: 0,
  };

  const coordinator: UiTransactionCoordinator = {
    events: state.events.asObservable(),
    begin() {
      if (state.disposed) {
        throw new Error('UI transaction coordinator is disposed');
      }
      state.sequence += 1;
      const transactionState: TransactionState = {
        coordinator: state,
        id: uiTransactionId(`tx-${String(state.sequence)}`),
        active: true,
        acknowledged: false,
        attachments: 0,
        releasedUnabsorbed: false,
      };
      const transaction: UiTransaction = {
        id: transactionState.id,
        acknowledge() {
          if (!transactionState.active || transactionState.acknowledged) {
            return;
          }
          transactionState.acknowledged = true;
          state.events.next({ type: 'acknowledged', id: transactionState.id });
        },
        reject() {
          finishTransaction(transactionState, 'rejected');
        },
        reconcile() {
          finishTransaction(transactionState, 'reconciled');
        },
      };
      transactionStates.set(transaction, transactionState);
      state.transactions.add(transactionState);
      return transaction;
    },
    reset() {
      if (state.disposed) {
        return;
      }
      for (const transaction of state.transactions) {
        transaction.active = false;
        transaction.attachments = 0;
      }
      state.transactions.clear();
      state.events.next({ type: 'reset' });
    },
    dispose() {
      if (state.disposed) {
        return;
      }
      state.disposed = true;
      for (const transaction of state.transactions) {
        transaction.active = false;
        transaction.attachments = 0;
      }
      state.transactions.clear();
      state.events.next({ type: 'disposed' });
      state.events.complete();
    },
  };
  coordinatorStates.set(coordinator, state);
  return coordinator;
}

export function createOptimisticProjection<T>(
  initialBase: T,
  coordinator: UiTransactionCoordinator,
): OptimisticProjection<T> {
  if (coordinatorStates.get(coordinator)?.disposed !== false) {
    throw new Error('A supported, active UI transaction coordinator is required');
  }

  let authoritativeBase = initialBase;
  // In order. An append pushes in place (O(1)); a removal makes a new array (`sweep`). It never escapes.
  let contributions: TrackedContribution<T>[] = [];
  /** Contributions held per transaction: a retirement of a transaction with none here is O(1). */
  const held = new Map<UiTransactionId, number>();
  /** The fold of `authoritativeBase` and every contribution, as last published: an append folds onto it. */
  let folded = initialBase;
  /** The contributions that could not apply in that fold, in order (a new array per change: reporting reads one). */
  let failing: readonly Failure<T>[] = [];
  let disposed = false;
  const projected = observableState(initialBase);
  const divergences = new Subject<OptimisticProjectionDivergence<T>>();

  const hold = (id: UiTransactionId, by: number) => {
    const count = (held.get(id) ?? 0) + by;
    if (count > 0) held.set(id, count);
    else held.delete(id);
  };

  // Publishes the fold, then reports each contribution that could not apply (still held: a handler may resolve one,
  // reject, and a contribution dropped meanwhile is not reported), so subscribers never see a partial fold.
  function publish(base: T): void {
    const report = failing;
    projected.set(folded);
    for (const { contribution, cause } of report) {
      if (contribution.dropped) continue;
      divergences.next({
        transactionId: contribution.transactionId,
        phase: 'application',
        authoritativeBase: base,
        cause,
      });
    }
  }

  // The base followed by every active contribution; one that cannot apply is skipped (and kept) instead of hiding the
  // ones after it.
  function emit(): void {
    const base = authoritativeBase;
    let value = base;
    const failed: Failure<T>[] = [];
    // As the list stands now: one appended meanwhile (re-entrantly) folds itself.
    const list = contributions;
    for (let i = 0, n = list.length; i < n; i++) {
      const contribution = list[i];
      try {
        value = contribution.apply(value);
      } catch (cause) {
        failed.push({ contribution, cause });
      }
    }
    folded = value;
    failing = failed;
    publish(base);
  }

  // One more contribution, last in order: folded onto the fold as it stands (the same value a refold gives, the
  // contributions before it unchanged), and the same reports (every contribution still failing, then this one). A
  // refold per contribution was K² for K in flight (view-models.md §5).
  function append(contribution: TrackedContribution<T>): void {
    contributions.push(contribution);
    hold(contribution.transactionId, 1);
    try {
      folded = contribution.apply(folded);
    } catch (cause) {
      failing = [...failing, { contribution, cause }];
    }
    publish(authoritativeBase);
  }

  /** Takes out the contributions marked `dropped`, in one pass. */
  function sweep(): void {
    contributions = contributions.filter((item) => !item.dropped);
  }

  function drop(transactionId: UiTransactionId): void {
    if (!held.has(transactionId)) {
      return;
    }
    held.delete(transactionId);
    const removed: TrackedContribution<T>[] = [];
    for (const item of contributions) {
      if (item.transactionId !== transactionId) continue;
      item.dropped = true;
      removed.push(item);
    }
    sweep();
    for (const item of removed) {
      item.release(false);
    }
    emit();
  }

  function clear(): void {
    const removed = contributions;
    contributions = [];
    held.clear();
    for (const item of removed) {
      item.dropped = true;
    }
    for (const item of removed) {
      item.release(false);
    }
    emit();
  }

  const eventSubscription = coordinator.events.subscribe((event) => {
    if (event.type === 'rejected' || event.type === 'reconciled') {
      drop(event.id);
      return;
    }
    if (event.type === 'reset') {
      clear();
      return;
    }
    if (event.type === 'disposed') {
      clear();
      disposed = true;
      divergences.complete();
    }
  });

  return {
    projected: projected.asReadonly,
    divergences: divergences.asObservable(),
    setBase(nextBase) {
      if (disposed) {
        return;
      }
      authoritativeBase = nextBase;
      // Decide first, then change the list as it stands, then report: a
      // divergence handler that rejects or contributes no longer has its change overwritten.
      const absorbed = new Set<TrackedContribution<T>>();
      const failed: Failure<T>[] = [];
      const list = contributions;
      for (let i = 0, n = list.length; i < n; i++) {
        const contribution = list[i];
        try {
          if (contribution.isAbsorbed?.(nextBase) === true) absorbed.add(contribution);
        } catch (cause) {
          failed.push({ contribution, cause });
        }
      }
      for (const contribution of absorbed) {
        contribution.dropped = true;
        hold(contribution.transactionId, -1);
      }
      if (absorbed.size > 0) sweep();
      emit();
      for (const contribution of absorbed) {
        contribution.release(true);
      }
      for (const { contribution, cause } of failed) {
        if (contribution.dropped) continue;
        divergences.next({
          transactionId: contribution.transactionId,
          phase: 'absorption',
          authoritativeBase: nextBase,
          cause,
        });
      }
    },
    contribute(transaction, contribution) {
      if (disposed) {
        throw new Error('Optimistic projection is disposed');
      }
      const attachment = attachTransaction(coordinator, transaction);
      append({
        transactionId: attachment.id,
        apply: contribution.apply,
        isAbsorbed: contribution.isAbsorbed,
        release: attachment.release,
        dropped: false,
      });
    },
    reset() {
      if (disposed) {
        return;
      }
      clear();
    },
    dispose() {
      if (disposed) {
        return;
      }
      clear();
      disposed = true;
      eventSubscription.unsubscribe();
      divergences.complete();
    },
  };
}
