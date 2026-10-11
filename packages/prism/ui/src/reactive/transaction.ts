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
};

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
  let contributions: readonly TrackedContribution<T>[] = [];
  let disposed = false;
  const projected = observableState(initialBase);
  const divergences = new Subject<OptimisticProjectionDivergence<T>>();

  // The base followed by every active contribution; one that cannot apply is
  // skipped (and kept) instead of hiding the ones after it, and divergences go out after the value,
  // so subscribers never see a partial fold. A handler may resolve one (reject): a contribution
  // dropped meanwhile is not reported.
  function emit(): void {
    const base = authoritativeBase;
    let value = base;
    const failed: { contribution: TrackedContribution<T>; cause: unknown }[] = [];
    for (const contribution of contributions) {
      try {
        value = contribution.apply(value);
      } catch (cause) {
        failed.push({ contribution, cause });
      }
    }
    projected.set(value);
    for (const { contribution, cause } of failed) {
      if (!contributions.includes(contribution)) continue;
      divergences.next({
        transactionId: contribution.transactionId,
        phase: 'application',
        authoritativeBase: base,
        cause,
      });
    }
  }

  function drop(transactionId: UiTransactionId): void {
    const removed = contributions.filter((item) => item.transactionId === transactionId);
    if (removed.length === 0) {
      return;
    }
    contributions = contributions.filter((item) => item.transactionId !== transactionId);
    for (const item of removed) {
      item.release(false);
    }
    emit();
  }

  function clear(): void {
    const removed = contributions;
    contributions = [];
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
      const failed: { contribution: TrackedContribution<T>; cause: unknown }[] = [];
      for (const contribution of contributions) {
        try {
          if (contribution.isAbsorbed?.(nextBase) === true) absorbed.add(contribution);
        } catch (cause) {
          failed.push({ contribution, cause });
        }
      }
      contributions = contributions.filter((contribution) => !absorbed.has(contribution));
      emit();
      for (const contribution of absorbed) {
        contribution.release(true);
      }
      for (const { contribution, cause } of failed) {
        if (!contributions.includes(contribution)) continue;
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
      contributions = [
        ...contributions,
        {
          transactionId: attachment.id,
          apply: contribution.apply,
          isAbsorbed: contribution.isAbsorbed,
          release: attachment.release,
        },
      ];
      emit();
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
