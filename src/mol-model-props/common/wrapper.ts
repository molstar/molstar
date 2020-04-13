/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifWriter } from '../../mol-io/writer/cif';
import { Model } from '../../mol-model/structure';
import { dateToUtcString } from '../../mol-util/date';
import { MmcifFormat } from '../../mol-model-formats/structure/mmcif';

interface PropertyWrapper<Data> {
    info: PropertyWrapper.Info,
    data: Data
}

namespace PropertyWrapper {
    export interface Info {
        timestamp_utc: string
    }

    export function createInfo(): Info {
        return { timestamp_utc: dateToUtcString(new Date()) };
    }

    export function defaultInfoCategory<Ctx>(name: string, getter: (ctx: Ctx) => Info | undefined): CifWriter.Category<Ctx> {
        return {
            name,
            instance(ctx) {
                const info = getter(ctx);
                return {
                    fields: _info_fields,
                    source: [{ data: info, rowCount: 1 }]
                };
            }
        };
    }

    const _info_fields: CifWriter.Field<number, Info>[] = [
        CifWriter.Field.str('updated_datetime_utc', (_, date) => date.timestamp_utc)
    ];

    export function tryGetInfoFromCif(categoryName: string, model: Model): Info | undefined {
        if (!MmcifFormat.is(model.sourceData) || !model.sourceData.data.frame.categoryNames.includes(categoryName)) {
            return;
        }

        const timestampField = model.sourceData.data.frame.categories[categoryName].getField('updated_datetime_utc');
        if (!timestampField || timestampField.rowCount === 0) return;

        return { timestamp_utc: timestampField.str(0) || dateToUtcString(new Date()) };
    }
}

export { PropertyWrapper };