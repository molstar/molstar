/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export interface PickingId {
  objectId: number;
  instanceId: number;
  groupId: number;
}

export namespace PickingId {
  /**
   * The null id of each target, its clear (white) read as 32 bits, 2^32 - 2: no object, no
   * instance; a group's: the whole instance (e.g. a clipped mesh's cap)
   */
  export const Null = 4294967294 as const;

  export function areSame(a: PickingId, b: PickingId) {
    return a.objectId === b.objectId && a.instanceId === b.instanceId && a.groupId === b.groupId;
  }
}
