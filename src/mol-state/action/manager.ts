/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateAction } from 'mol-state/action';
import { StateObject } from '../object';

export { StateActionManager }

class StateActionManager {
    private actions: Map<StateAction['id'], StateAction> = new Map();
    private fromTypeIndex = new Map<StateObject.Type, StateAction[]>();

    add(action: StateAction) {
        if (this.actions.has(action.id)) return this;

        this.actions.set(action.id, action);

        for (const t of action.definition.from) {
            if (this.fromTypeIndex.has(t.type)) {
                this.fromTypeIndex.get(t.type)!.push(action);
            } else {
                this.fromTypeIndex.set(t.type, [action]);
            }
        }

        return this;
    }

    fromType(type: StateObject.Type): ReadonlyArray<StateAction> {
        return this.fromTypeIndex.get(type) || [];
    }
}