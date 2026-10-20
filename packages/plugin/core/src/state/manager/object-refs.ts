/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { State, StateObject } from '@molstar/core/state';
import type { PluginComponent } from '../component.js';

/**
 * The refs of a state's cells whose object passes `test`, kept from the object events (seeded from the cells). A ref
 * may stay after its object was cleared within a pass (created and removed in one update): a superset, never fewer.
 * A hierarchy manager that tracks none of its root objects skips walking the whole state on every change.
 */
export function trackObjectRefs(
  component: PluginComponent,
  state: State,
  test: (obj: StateObject | undefined) => boolean,
): ReadonlySet<string> {
  const refs = new Set<string>();
  state.cells.forEach((c, ref) => {
    if (test(c.obj)) refs.add(ref);
  });
  component.subscribe(state.events.object.created, (e) => {
    if (test(e.obj)) refs.add(e.ref);
  });
  component.subscribe(state.events.object.updated, (e) => {
    if (test(e.obj)) refs.add(e.ref);
    else refs.delete(e.ref);
  });
  component.subscribe(state.events.object.removed, (e) => refs.delete(e.ref));
  return refs;
}
