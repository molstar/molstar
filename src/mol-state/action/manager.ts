/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateAction } from '../action';
import { StateObject, StateObjectCell } from '../object';
import { StateTransformer } from '../transformer';

export { StateActionManager }

class StateActionManager {
    private actions: Map<StateAction['id'], StateAction> = new Map();
    private fromTypeIndex = new Map<StateObject.Type, StateAction[]>();

    add(actionOrTransformer: StateAction | StateTransformer) {
        const action = StateTransformer.is(actionOrTransformer) ? actionOrTransformer.toAction() : actionOrTransformer;

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

    fromCell(cell: StateObjectCell, ctx: unknown): ReadonlyArray<StateAction> {
        const obj = cell.obj;
        if (!obj) return [];

        const actions = this.fromTypeIndex.get(obj.type);
        if (!actions) return [];
        let hasTest = false;
        for (const a of actions) {
            if (a.definition.isApplicable) {
                hasTest = true;
                break;
            }
        }
        if (!hasTest) return actions;

        const ret: StateAction[] = [];
        for (const a of actions) {
            if (a.definition.isApplicable) {
                if (a.definition.isApplicable(obj, cell.transform, ctx)) {
                    ret.push(a);
                }
            } else {
                ret.push(a);
            }
        }
        return ret;
    }
}