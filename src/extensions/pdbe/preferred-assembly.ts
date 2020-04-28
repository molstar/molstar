/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, Table } from '../../mol-data/db';
import { toTable } from '../../mol-io/reader/cif/schema';
import { CifWriter } from '../../mol-io/writer/cif';
import { Model } from '../../mol-model/structure';
import { ModelSymmetry } from '../../mol-model-formats/structure/property/symmetry';
import { MmcifFormat } from '../../mol-model-formats/structure/mmcif';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';

export namespace PDBePreferredAssembly {
    export type Property = string

    export function getFirstFromModel(model: Model): Property {
        const symmetry = ModelSymmetry.Provider.get(model);
        return symmetry?.assemblies.length ? symmetry.assemblies[0].id : '';
    }

    export function get(model: Model): Property {
        return model._staticPropertyData.__PDBePreferredAssebly__ || getFirstFromModel(model);
    }
    function set(model: Model, prop: Property) {
        (model._staticPropertyData.__PDBePreferredAssebly__ as Property) = prop;
    }

    export const Schema = {
        pdbe_preferred_assembly: {
            assembly_id: Column.Schema.str
        }
    };
    export type Schema = typeof Schema

    export const Descriptor = CustomPropertyDescriptor({
        name: 'pdbe_preferred_assembly',
        cifExport: {
            prefix: 'pdbe',
            context(ctx): Property { return get(ctx.firstModel); },
            categories: [{
                name: 'pdbe_preferred_assembly',
                instance(ctx: Property) {
                    return CifWriter.Category.ofTable(Table.ofArrays(Schema.pdbe_preferred_assembly, { assembly_id: [ctx] }));
                }
            }]
        }
    });

    function fromCifData(model: Model): string | undefined {
        if (!MmcifFormat.is(model.sourceData)) return void 0;
        const cat = model.sourceData.data.frame.categories.pdbe_preferred_assembly;
        if (!cat) return void 0;
        return toTable(Schema.pdbe_preferred_assembly, cat).assembly_id.value(0) || getFirstFromModel(model);
    }

    export async function attachFromCifOrApi(model: Model, params: {
        // optional JSON source
        PDBe_apiSourceJson?: (model: Model) => Promise<any>
    }) {
        if (model.customProperties.has(Descriptor)) return true;

        let asmName: string | undefined = fromCifData(model);
        if (asmName === void 0 &&  params.PDBe_apiSourceJson) {
            const data = await params.PDBe_apiSourceJson(model);
            if (!data) return false;
            asmName = asmNameFromJson(model, data);
        } else {
            return false;
        }

        if (!asmName) return false;

        model.customProperties.add(Descriptor);
        set(model, asmName);
        return true;
    }
}

function asmNameFromJson(modelData: Model, data: any): string {
    const assemblies = data[0] && data[0].assemblies;
    if (!assemblies || !assemblies.length) return PDBePreferredAssembly.getFirstFromModel(modelData);

    for (const asm of assemblies) {
        if (asm.preferred) {
            return asm.assembly_id;
        }
    }
    return assemblies[0].assembly_id;
}