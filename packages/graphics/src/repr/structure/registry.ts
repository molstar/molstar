/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Structure } from '@molstar/model/model/structure';
import { objectForEach } from '@molstar/core/util/object';
import { RepresentationRegistry, type RepresentationProvider } from '../representation.js';
import { BuiltInStructureRepresentations } from './catalog.js';
import type { StructureRepresentationState } from './representation.js';

export class StructureRepresentationRegistry extends RepresentationRegistry<Structure, StructureRepresentationState> {
  constructor() {
    super();
    objectForEach(BuiltInStructureRepresentations, (p) => this.add(p as any));
  }
}

export namespace StructureRepresentationRegistry {
  type _BuiltIn = typeof BuiltInStructureRepresentations;
  export type BuiltIn = keyof _BuiltIn;
  export type BuiltInParams<T extends BuiltIn> = Partial<RepresentationProvider.ParamValues<_BuiltIn[T]>>;
}
