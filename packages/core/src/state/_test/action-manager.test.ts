/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateAction, StateObject, type StateObjectCell, StateTransformer } from '../index.js';
import { StateActionManager } from '../action/manager.js';

interface TypeInfo {
  name: string;
  typeClass: 'Root' | 'Data';
}
const Create = StateObject.factory<TypeInfo>();

class Root extends Create({ name: 'Root', typeClass: 'Root' }) {}
class Other extends Create({ name: 'Other', typeClass: 'Data' }) {}

const NS = 'action-manager-spec';
let counter = 0;

function makeTransformer() {
  return StateTransformer.create<Root, Other, {}>(NS, {
    name: `to-other-${counter++}`,
    from: [Root],
    to: [Other],
    display: { name: 'To Other' },
    apply() {
      return new Other({ name: 'x', typeClass: 'Data' });
    },
  });
}

function makeAction() {
  return StateAction.create<Root, void, {}>({
    from: [Root],
    display: { name: 'Action' },
    run() {},
  });
}

function rootCell() {
  return { obj: new Root({ name: 'Root', typeClass: 'Root' }), transform: {} } as unknown as StateObjectCell;
}

describe('StateTransformer.toAction', () => {
  it('returns one cached action per transformer', () => {
    const t = makeTransformer();
    const a = t.toAction();
    expect(t.toAction()).toBe(a);
    expect(t.toAction().id).toBe(a.id);
    expect(makeTransformer().toAction()).not.toBe(a);
  });
});

describe('StateActionManager', () => {
  function setup() {
    const m = new StateActionManager();
    const added = jest.fn();
    const removed = jest.fn();
    m.events.added.subscribe(added);
    m.events.removed.subscribe(removed);
    return { m, added, removed };
  }

  it('counts registrations and emits only when presence changes', () => {
    const { m, added, removed } = setup();
    const action = makeAction();
    const cell = rootCell();

    m.add(action);
    m.add(action);
    expect(added).toHaveBeenCalledTimes(1);
    expect(m.fromCell(cell, undefined)).toEqual([action]);

    m.remove(action);
    expect(removed).toHaveBeenCalledTimes(0);
    expect(m.fromCell(cell, undefined)).toEqual([action]);

    m.remove(action);
    expect(removed).toHaveBeenCalledTimes(1);
    expect(m.fromCell(cell, undefined)).toEqual([]);

    m.add(action);
    expect(added).toHaveBeenCalledTimes(2);
    expect(m.fromCell(cell, undefined)).toEqual([action]);
  });

  it('treats a transformer and its action as the same registration', () => {
    const { m, added, removed } = setup();
    const t = makeTransformer();
    const cell = rootCell();

    m.add(t);
    m.add(t.toAction());
    m.add(t);
    expect(added).toHaveBeenCalledTimes(1);
    expect(m.fromCell(cell, undefined)).toEqual([t.toAction()]);

    m.remove(t.toAction().id);
    m.remove(t);
    expect(removed).toHaveBeenCalledTimes(0);
    m.remove(t.toAction());
    expect(removed).toHaveBeenCalledTimes(1);
    expect(m.fromCell(cell, undefined)).toEqual([]);
  });

  it('treats unknown removal as a no-op', () => {
    const { m, added, removed } = setup();
    const action = makeAction();
    const other = makeAction();

    m.remove(action);
    m.remove(makeTransformer());
    m.remove(action.id);
    expect(removed).toHaveBeenCalledTimes(0);

    m.add(action);
    m.remove(other);
    expect(removed).toHaveBeenCalledTimes(0);
    expect(m.fromCell(rootCell(), undefined)).toEqual([action]);

    m.remove(action);
    m.remove(action);
    expect(removed).toHaveBeenCalledTimes(1);

    // extra removals did not leave a negative count
    m.add(action);
    expect(added).toHaveBeenCalledTimes(2);
    m.remove(action);
    expect(removed).toHaveBeenCalledTimes(2);
  });

  it('never reports a conflict', () => {
    const { m } = setup();
    const action = makeAction();
    const t = makeTransformer();
    expect(m.findConflict(action)).toBeUndefined();
    m.add(action);
    m.add(t);
    expect(m.findConflict(action)).toBeUndefined();
    expect(m.findConflict(t)).toBeUndefined();
  });
});
