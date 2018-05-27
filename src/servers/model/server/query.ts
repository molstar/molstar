/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from 'mol-util';

export interface ResponseFormat {

}

export interface Query {
    id: UUID,

    sourceId: 'file' | string,
    entryId: string,

    kind: string,
    params: any,

    responseFormat: ResponseFormat
}

// export class QueryQueue {

// }

export function resolveQuery() {

}