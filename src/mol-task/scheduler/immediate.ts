/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

namespace ImmediateScheduler {
    // Adds the function to the start of the "immediate queue"
    export async function first<T>(f: () => T): Promise<T> {
        return f();
    }

    // Adds the function to the end of the "immediate queue"
    export async function last<T>(f: () => T): Promise<T> {
        return f();
    }
}

export default ImmediateScheduler