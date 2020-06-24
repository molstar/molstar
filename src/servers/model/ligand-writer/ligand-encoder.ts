/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { StringBuilder } from '../../../mol-util';
import Writer from '../../../mol-io/writer/writer';
import { Encoder, Category } from '../../../mol-io/writer/cif/encoder';
import { ComponentBond } from '../../../mol-model-formats/structure/property/bonds/comp';

export abstract class LigandExplorer implements Encoder<string> {
    protected builder: StringBuilder;
    protected componentData: ComponentBond;
    readonly isBinary = false;
    binaryEncodingProvider = void 0;

    abstract writeCategory<Ctx>(category: Category<Ctx>, context?: Ctx): void;

    abstract encode(): void;

    setComponentBondData(componentData: ComponentBond) {
        this.componentData = componentData;
    }

    writeTo(stream: Writer) {
        const chunks = StringBuilder.getChunks(this.builder);
        for (let i = 0, _i = chunks.length; i < _i; i++) {
            stream.writeString(chunks[i]);
        }
    }

    getSize() {
        return StringBuilder.getSize(this.builder);
    }

    getData() {
        return StringBuilder.getString(this.builder);
    }

    protected skipHydrogen(label: string) {
        if (this.hydrogens) {
            return false;
        }
        return label.startsWith('H');
    }

    protected getLabel(s: string) {
        return s.replace(/[^A-Z]+/g, '');
    }

    protected getSortedFields<Ctx>(instance: Category.Instance<Ctx>, names: string[]) {
        return names.map(n => this.getField(instance, n));
    }

    protected getField<Ctx>(instance: Category.Instance<Ctx>, name: string) {
        return instance.fields.find(f => f.name === name)!;
    }

    startDataBlock() {}

    setFilter() {}

    setFormatter() {}

    isCategoryIncluded() {
        return true;
    }

    constructor(readonly encoder: string, readonly hydrogens: boolean) {
        this.builder = StringBuilder.create();
    }
}