/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Annotation } from './annotation/annotation'
import { UUID } from 'mol-util'

interface Annotations {
    id: UUID,
    all: Annotation[],
    byKind: { [kind: number]: Annotation }
    //getAll()
}

export { Annotations }