/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** A mutable value cell. */
interface ValueCell<T> { value: T }
/** Create a mutable value cell. */
function ValueCell<T>(value: T): ValueCell<T> { return { value }; }

/** An immutable value box that also holds a version of the attribute. */
interface ValueBox<T> { readonly version: number, readonly value: T }
/** Create a new box with the specified value and version = 0 */
function ValueBox<T>(value: T): ValueBox<T>
/** Create a new box by updating the value of an old box and incrementing the version number. */
function ValueBox<T>(box: ValueBox<T>, value: T): ValueBox<T>
function ValueBox<T>(boxOrValue: T | ValueBox<T>, value?: T): ValueBox<T> {
    if (arguments.length === 2) return { version: (boxOrValue as ValueBox<T>).version + 1, value: value! };
    return { version: 0, value: boxOrValue as T };
}

export { ValueCell, ValueBox };

