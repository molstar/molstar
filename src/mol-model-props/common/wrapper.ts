/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifWriter } from 'mol-io/writer/cif';
import { Model } from 'mol-model/structure';
import { CifExportContext } from 'mol-model/structure';
import { dateToUtcString } from 'mol-util/date';

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

    export function defaultInfoCategory(name: string, getter: (model: Model) => PropertyWrapper<unknown> | undefined): CifWriter.Category<CifExportContext> {
        return {
            name,
            instance(ctx) {
                const prop = getter(ctx.firstModel);
                if (!prop) return CifWriter.Category.Empty;
                return {
                    fields: _info_fields,
                    source: [{ data: prop.info.timestamp_utc, rowCount: 1 }]
                };
            }
        }
    }

    const _info_fields: CifWriter.Field<number, string>[] = [
        CifWriter.Field.str('updated_datetime_utc', (_, date) => date)
    ];

    export function tryGetInfoFromCif(categoryName: string, model: Model): Info | undefined {
        if (model.sourceData.kind !== 'mmCIF' || !model.sourceData.frame.categoryNames.includes(categoryName)) {
            return;
        }

        const timestampField = model.sourceData.frame.categories[categoryName].getField('updated_datetime_utc');
        if (!timestampField || timestampField.rowCount === 0) return;

        return { timestamp_utc: timestampField.str(0) || dateToUtcString(new Date()) };
    }
}

export { PropertyWrapper }