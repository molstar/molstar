/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StateAction } from '../action.js';
import type { StateObject, StateObjectCell } from '../object.js';
import { StateTransformer } from '../transformer.js';
import { UUID } from '@molstar/core/util';
import { arraySetRemove } from '@molstar/core/util/array';
import { RxEventHelper } from '@molstar/core/util/rx-event-helper';

export { StateActionManager };

class StateActionManager {
  private ev = RxEventHelper.create();
  /** Registered actions with the number of registrations, keyed by action id. */
  private actions: Map<StateAction['id'], { action: StateAction; count: number }> = new Map();
  private fromTypeIndex = new Map<StateObject.Type, StateAction[]>();

  readonly events = {
    added: this.ev<undefined>(),
    removed: this.ev<undefined>(),
  };

  /** Registering the same action again increments its count; `added` fires only when the action first appears. */
  add(actionOrTransformer: StateAction | StateTransformer) {
    const action = StateTransformer.is(actionOrTransformer) ? actionOrTransformer.toAction() : actionOrTransformer;

    const entry = this.actions.get(action.id);
    if (entry) {
      entry.count++;
      return this;
    }

    this.actions.set(action.id, { action, count: 1 });

    for (const t of action.definition.from) {
      if (this.fromTypeIndex.has(t.type)) {
        this.fromTypeIndex.get(t.type)!.push(action);
      } else {
        this.fromTypeIndex.set(t.type, [action]);
      }
    }

    this.events.added.next(void 0);

    return this;
  }

  /**
   * Decrements the count of the action; it is removed (and `removed` fires) when the count reaches zero.
   * Removing an unknown action is a no-op.
   */
  remove(actionOrTransformer: StateAction | StateTransformer | UUID) {
    const id = StateTransformer.is(actionOrTransformer)
      ? actionOrTransformer.toAction().id
      : UUID.is(actionOrTransformer)
        ? actionOrTransformer
        : actionOrTransformer.id;

    const entry = this.actions.get(id);
    if (!entry) return this;

    if (--entry.count > 0) return this;

    const { action } = entry;
    this.actions.delete(id);
    for (const t of action.definition.from) {
      const xs = this.fromTypeIndex.get(t.type);
      if (!xs) continue;

      arraySetRemove(xs, action);
      if (xs.length === 0) this.fromTypeIndex.delete(t.type);
    }

    this.events.removed.next(void 0);

    return this;
  }

  /** Registering an action never conflicts: the same action is counted, distinct actions have distinct ids. */
  findConflict(_actionOrTransformer: StateAction | StateTransformer): string | undefined {
    return undefined;
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

  dispose() {
    this.ev.dispose();
  }
}
