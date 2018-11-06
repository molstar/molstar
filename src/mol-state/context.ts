/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject } from './object';
import { Transform } from './transform';
import { RxEventHelper } from 'mol-util/rx-event-helper';

export { StateContext }

class StateContext {
    private ev = RxEventHelper.create();

    readonly events = {
        object: {
            stateChanged: this.ev<{ ref: Transform.Ref }>(),
            propsChanged: this.ev<{ ref: Transform.Ref, newProps: unknown }>(),

            updated: this.ev<{ ref: Transform.Ref, obj?: StateObject }>(),
            replaced: this.ev<{ ref: Transform.Ref, oldObj?: StateObject, newObj?: StateObject }>(),
            created: this.ev<{ ref: Transform.Ref, obj: StateObject }>(),
            removed: this.ev<{ ref: Transform.Ref, obj?: StateObject }>(),

            currentChanged: this.ev<{ ref: Transform.Ref }>()
        },
        warn: this.ev<string>(),
        updated: this.ev<void>()
    };

    readonly behaviours = {
        currentObject: this.ev.behavior<{ ref: Transform.Ref }>(void 0 as any)
    };

    readonly globalContext: unknown;
    readonly defaultObjectProps: unknown;

    dispose() {
        this.ev.dispose();
    }

    constructor(params: { globalContext: unknown, defaultObjectProps: unknown, rootRef: Transform.Ref }) {
        this.globalContext = params.globalContext;
        this.defaultObjectProps = params.defaultObjectProps;
        this.behaviours.currentObject.next({ ref: params.rootRef });
    }
}