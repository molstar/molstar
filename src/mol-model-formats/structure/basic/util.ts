/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BasicData } from './schema';
import { Table } from '../../../mol-data/db';

export function getModelGroupName(model_id: number, data: BasicData) {
    const { ihm_model_group, ihm_model_group_link } = data;

    const link = Table.pickRow(ihm_model_group_link, i => ihm_model_group_link.model_id.value(i) === model_id);
    if (link) {
        const group = Table.pickRow(ihm_model_group, i => ihm_model_group.id.value(i) === link.group_id);
        if (group) return group.name;
    }
    return '';
}