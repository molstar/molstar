/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** A mutable value reference. */
interface ValueRef<T> { ref: T }

namespace ValueRef {
    export function create<T>(ref: T): ValueRef<T> { return { ref }; }
    export function set<T>(ref: ValueRef<T>, value: T) { ref.ref = value; return ref; }
}

let _valueBoxId = 0;
function getNextId() {
    return _valueBoxId++ % 0x7FFFFFFF;
}

/**
 * An immutable value box that also holds a version of the attribute.
 * Optionally includes automatically propadated "metadata".
 */
type ValueBox<T, D = never> = {
    // Unique identifier in the range 0 to 0x7FFFFFFF
    readonly id: number,
    readonly version: number,
    readonly metadata: D,
    readonly value: T,
}

namespace ValueBox {
    export function create<T, D = never>(value: T, metadata?: D): ValueBox<T, D> {
        return { id: getNextId(), version: 0, value, metadata: metadata! };
    }

    /** If diffInfo is not specified, copy the old value */
    export function withValue<T, D>(box: ValueBox<T, D>, value: T): ValueBox<T, D> {
        return { id: box.id, version: box.version + 1, value, metadata: box.metadata };
    }
}

/** An immutable box stored inside a mutable cell. */
type ValueCell<T, D = never> = ValueRef<ValueBox<T, D>>

namespace ValueCell {
    export function create<T, D = never>(value: T, metadata?: D): ValueCell<T, D> {
        return ValueRef.create(ValueBox.create(value, metadata));
    }

    /** If diffInfo is not specified, copy the old value */
    export function update<T, D>(cell: ValueCell<T, D>, value: T): ValueCell<T, D> {
        return ValueRef.set(cell, ValueBox.withValue(cell.ref, value));
    }

    export function set<T, D>(cell: ValueCell<T, D>, box: ValueBox<T, D>): ValueCell<T, D> {
        return ValueRef.set(cell, box);
    }
}

export { ValueRef, ValueBox, ValueCell };