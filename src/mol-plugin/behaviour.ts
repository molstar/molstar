/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export { PluginBehaviour }

interface PluginBehaviour<P> {
    register(): void,
    unregister(): void,

    /** Update params in place. Optionally return a promise if it depends on an async action. */
    update(params: P): void | Promise<void>
}

namespace PluginBehaviour {
    export interface Ctor<P> {
        create(params: P): PluginBehaviour<P>
    }
}