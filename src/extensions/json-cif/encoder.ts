/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../mol-data/db';
import { Category, Encoder } from '../../mol-io/writer/cif/encoder';
import { BinaryEncodingProvider } from '../../mol-io/writer/cif/encoder/binary';
import { getCategoryInstanceData, getIncludedFields } from '../../mol-io/writer/cif/encoder/util';
import { Writer } from '../../mol-io/writer/writer';
import { JSONCifCategory, JSONCifDataBlock, JSONCifFile, JSONCifVERSION } from './model';

export class JSONCifEncoder implements Encoder<string> {
    private data: JSONCifFile;
    private dataBlocks: JSONCifDataBlock[] = [];
    private encodedData: string | undefined;
    private filter: Category.Filter = Category.DefaultFilter;

    readonly isBinary = false;
    readonly binaryEncodingProvider: BinaryEncodingProvider | undefined;

    setFilter(filter?: Category.Filter) {
        this.filter = filter || Category.DefaultFilter;
    }

    isCategoryIncluded(name: string) {
        return this.filter.includeCategory(name);
    }

    setFormatter(formatter?: Category.Formatter) {
        // No formatter needed for JSON encoding.
    }

    startDataBlock(header: string) {
        this.dataBlocks.push({
            header: (header || '').replace(/[ \n\t]/g, '').toUpperCase(),
            categoryNames: [],
            categories: {}
        });
    }

    writeCategory<Ctx>(category: Category<Ctx>, context?: Ctx, options?: Encoder.WriteCategoryOptions) {
        if (this.encodedData) {
            throw new Error('The writer contents have already been encoded, no more writing.');
        }

        if (!this.dataBlocks.length) {
            throw new Error('No data block created.');
        }

        if (!options?.ignoreFilter && !this.filter.includeCategory(category.name)) return;

        const { instance, rowCount, source } = getCategoryInstanceData(category, context);
        if (!rowCount) return;

        const fields = getIncludedFields(instance);
        if (!fields.length) return;

        const rows: Record<string, any>[] = [];
        const cat: JSONCifCategory = { name: category.name, fieldNames: fields.map(f => f.name), rows };

        for (const src of source) {
            const d = src.data;
            const keys = src.keys();
            while (keys.hasNext) {
                const row: Record<string, any> = {};
                const k = keys.move();
                for (const f of fields) {
                    const kind = f.valueKind ? f.valueKind(k, d) : Column.ValueKinds.Present;
                    if (kind === Column.ValueKinds.Present) {
                        row[f.name] = f.value(k, d, rows.length);
                    } else if (kind === Column.ValueKinds.NotPresent) {
                        row[f.name] = null;
                    }
                }
                cat.rows.push(row);
            }
        }

        this.dataBlocks[this.dataBlocks.length - 1].categoryNames.push(cat.name);
        this.dataBlocks[this.dataBlocks.length - 1].categories[cat.name] = cat;
    }

    encode() {
        if (this.encodedData) return;
        this.encodedData = this.options?.formatJSON ? JSON.stringify(this.data, null, 2) : JSON.stringify(this.data);
    }

    writeTo(writer: Writer) {
        writer.writeString(this.encodedData!);
    }

    getData() {
        this.encode();
        return this.encodedData!;
    }

    getSize() {
        return this.encodedData?.length ?? 0;
    }

    getFile() {
        return this.data;
    }

    constructor(encoder: string, public options?: { formatJSON?: boolean }) {
        this.data = {
            encoder,
            version: JSONCifVERSION,
            dataBlocks: this.dataBlocks
        };
    }
}