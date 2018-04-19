/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, ElementSet, Element } from '../structure'

interface Annotation<E = any> {
    definition: Annotation.Definition<E>,
    getValue(l: Element.Location): E | undefined,
    getAll(l: ElementSet): { annotations: E[], /* TODO: map annotations to elements */ }
}

namespace Annotation {
    export const enum Kind {
        Atom,
        Residue,
        Sequence,
        Chain,
        Entity,
        Coarse,
        Spatial
    }

    export const enum Type {
        Num,
        Str,
        Obj
    }

    export interface Definition<E = any> {
        name: string,
        kind: Kind,
        type: Type,
        prepare<Data>(s: Structure, data: Data): Annotation<E>,
    }
}

export { Annotation }