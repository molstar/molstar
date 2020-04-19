/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { QueryName, QueryParams } from './api';

export interface MultipleQueryEntry<Name extends QueryName = QueryName> {
    data_source?: string,
    entryId: string,
    query: Name,
    params?: QueryParams<Name>,
    model_nums?: number[],
    copy_all_categories?: boolean
}

export interface MultipleQuerySpec {
    queries: MultipleQueryEntry[],
    encoding?: 'cif' | 'bcif',
    asTarGz?: boolean
}

export function getMultiQuerySpecFilename() {
    const date = new Date();
    return `result_${date.getMonth() + 1}-${date.getDate()}-${date.getHours()}-${date.getMinutes()}-${date.getSeconds()}.tar.gz`;
}