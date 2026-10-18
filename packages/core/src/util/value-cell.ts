/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from './id-factory.js';

/** A mutable value reference. */
interface ValueRef<T> {
  ref: T;
  /** Flagged when the reference changes through `ValueCell.update` or `ValueCell.set`, see `ValueCell.observe` */
  observers?: Set<ValueCell.Observer>;
  /**
   * For cells derived from other state (e.g. a value depending on the current render target): the value may change
   * without notifying the observers, which then check the version every time
   */
  readonly volatile?: boolean;
}

namespace ValueRef {
  export function create<T>(ref: T): ValueRef<T> {
    return { ref };
  }
  export function set<T>(ref: ValueRef<T>, value: T) {
    ref.ref = value;
    return ref;
  }
}

const getNextId = idFactory(0, 0x7fffffff);

/**
 * An immutable value box that also holds a version of the attribute.
 * Optionally includes automatically propadated "metadata".
 */
type ValueBox<T, D = never> = {
  /** Unique identifier in the range 0 to 0x7FFFFFFF */
  readonly id: number;
  readonly version: number;
  readonly metadata: D;
  readonly value: T;
};

namespace ValueBox {
  export function create<T, D = never>(value: T, metadata?: D): ValueBox<T, D> {
    return { id: getNextId(), version: 0, value, metadata: metadata! };
  }

  /** The box.metadata is carried over from the old box */
  export function withValue<T, D>(box: ValueBox<T, D>, value: T): ValueBox<T, D> {
    return { id: box.id, version: box.version + 1, value, metadata: box.metadata };
  }
}

/** An immutable box stored inside a mutable cell. */
type ValueCell<T, D = never> = ValueRef<ValueBox<T, D>>;

namespace ValueCell {
  export function create<T, D = never>(value: T, metadata?: D): ValueCell<T, D> {
    return ValueRef.create(ValueBox.create(value, metadata));
  }

  /**
   * Set to `changed` when an observed cell changes; lets the observer skip checking the versions of its cells while
   * none changed (e.g. the uniform block writers, per draw)
   */
  export interface Observer {
    changed: boolean;
  }

  /** `observer.changed` is set when `cell` changes, until `unobserve` */
  export function observe(cell: ValueCell<any, any>, observer: Observer) {
    (cell.observers ??= new Set()).add(observer);
  }

  export function unobserve(cell: ValueCell<any, any>, observer: Observer) {
    cell.observers?.delete(observer);
  }

  function notify(cell: ValueCell<any, any>) {
    const observers = cell.observers;
    if (observers) for (const o of observers) o.changed = true;
  }

  /** The box.metadata is carried over from the old box */
  export function update<T, D>(cell: ValueCell<T, D>, value: T): ValueCell<T, D> {
    ValueRef.set(cell, ValueBox.withValue(cell.ref, value));
    notify(cell);
    return cell;
  }

  export function set<T, D>(cell: ValueCell<T, D>, box: ValueBox<T, D>): ValueCell<T, D> {
    ValueRef.set(cell, box);
    notify(cell);
    return cell;
  }

  /** Updates the cell if the value is has changed, comparing by reference */
  export function updateIfChanged<T, D>(cell: ValueCell<T, D>, value: T): ValueCell<T, D> {
    return cell.ref.value !== value ? update(cell, value) : cell;
  }
}

export { ValueRef, ValueBox, ValueCell };
