/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, Table } from '../../mol-data/db';
import { toTable } from '../../mol-io/reader/cif/schema';
import { Model, CustomPropertyDescriptor } from '../../mol-model/structure';
import { mmCIF_chemCompBond_schema } from '../../mol-io/reader/cif/schema/mmcif-extras';
import { CifWriter } from '../../mol-io/writer/cif';

export namespace ChemCompBond {
    export type Property = Table<Schema['chem_comp_bond']>

    export function getFromModel(model: Model): Property {
        if (model.sourceData.kind !== 'mmCIF') return Table.ofUndefinedColumns(Schema.chem_comp_bond, 0);
        const { chem_comp_bond } = model.sourceData.data
        return Table.ofColumns(Schema.chem_comp_bond, {
            ...chem_comp_bond,
            molstar_protonation_variant: Column.Undefined(chem_comp_bond._rowCount, Column.Schema.Str())
        });
    }

    export function get(model: Model): Property {
        return model._staticPropertyData.__ChemCompBond__ || getFromModel(model);
    }
    function set(model: Model, prop: Property) {
        (model._staticPropertyData.__ChemCompBond__ as Property) = prop;
    }

    export const Schema = { chem_comp_bond: mmCIF_chemCompBond_schema };
    export type Schema = typeof Schema

    export const Descriptor = CustomPropertyDescriptor({
        isStatic: true,
        name: 'chem_comp_bond',
        cifExport: {
            prefix: '',
            context(ctx): Property { return get(ctx.firstModel); },
            categories: [{
                name: 'chem_comp_bond',
                instance(ctx: Property) {
                    return CifWriter.Category.ofTable(ctx);
                }
            }]
        }
    });

    function fromCifData(model: Model): Table<Schema['chem_comp_bond']> | undefined {
        if (model.sourceData.kind !== 'mmCIF') return void 0;
        const cat = model.sourceData.frame.categories.chem_comp_bond;
        if (!cat) return void 0;
        return toTable(Schema.chem_comp_bond, cat);
    }

    export async function attachFromCifOrTable(model: Model, params: {
        // optional Table source
        wwPDB_apiSourceTable?: (model: Model) => Promise<Table<Schema['chem_comp_bond']>>
    }) {
        if (model.customProperties.has(Descriptor)) return true;

        let chemCompBond: Table<Schema['chem_comp_bond']> | undefined = fromCifData(model);
        if (chemCompBond === void 0 && params.wwPDB_apiSourceTable) {
            const data = await params.wwPDB_apiSourceTable(model);
            if (!data) return false;
            chemCompBond = chemCompBondFromTable(model, data);
        } else {
            return false;
        }

        if (!chemCompBond) return false;

        model.customProperties.add(Descriptor);
        set(model, chemCompBond);
        return true;
    }
}

function chemCompBondFromTable(model: Model, table: Table<ChemCompBond.Schema['chem_comp_bond']>): Table<ChemCompBond.Schema['chem_comp_bond']> {
    return Table.pick(table, ChemCompBond.Schema.chem_comp_bond, (i: number) => {
        return model.properties.chemicalComponentMap.has(table.comp_id.value(i))
    })
}