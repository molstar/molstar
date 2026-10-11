/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, Table } from '@molstar/core/data/db';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { RuntimeContext } from '@molstar/core/task';
import { type GraphicsContext } from '@molstar/graphics/gl/context';
import { Representation, type RepresentationProvider } from '@molstar/graphics/repr/representation';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { Theme } from '@molstar/graphics/theme/theme';
import { createModels } from '@molstar/model/formats/structure/basic/parser';
import { BasicSchema, createBasic } from '@molstar/model/formats/structure/basic/schema';
import { ComponentBuilder } from '@molstar/model/formats/structure/common/component';
import { EntityBuilder } from '@molstar/model/formats/structure/common/entity';
import { Structure } from '@molstar/model/model/structure';
import { MoleculeType } from '@molstar/model/model/structure/model/types';

let themeRegistries: { colorThemeRegistry: ColorTheme.Registry; sizeThemeRegistry: SizeTheme.Registry } | undefined;
function getThemeRegistries() {
    if (!themeRegistries) {
        themeRegistries = {
            colorThemeRegistry: ColorTheme.createRegistry(),
            sizeThemeRegistry: SizeTheme.createRegistry(),
        };
    }
    return themeRegistries;
}

/** Two chains of carbon atoms on helices around the z axis, one atom per residue */
export async function helixStructure(atomsPerChain: number): Promise<Structure> {
    const chains = ['A', 'B'];
    const count = chains.length * atomsPerChain;
    const x = new Float32Array(count),
        y = new Float32Array(count),
        z = new Float32Array(count);
    const asymIds: string[] = [];
    const seqIds: number[] = [];
    for (let c = 0; c < chains.length; ++c) {
        for (let i = 0; i < atomsPerChain; ++i) {
            const k = c * atomsPerChain + i;
            const a = i * 0.6 + c * Math.PI;
            x[k] = 4 * Math.cos(a);
            y[k] = 4 * Math.sin(a);
            z[k] = i * 0.5 - atomsPerChain * 0.25;
            asymIds.push(chains[c]);
            seqIds.push(i + 1);
        }
    }

    const auth_asym_id = Column.ofStringArray(asymIds);
    const auth_atom_id = Column.ofConst('CA', count, Column.Schema.str);
    const auth_comp_id = Column.ofConst('ALA', count, Column.Schema.str);
    const auth_seq_id = Column.ofIntArray(seqIds);
    const atom_site = Table.ofPartialColumns(
        BasicSchema.atom_site,
        {
            auth_asym_id,
            auth_atom_id,
            auth_comp_id,
            auth_seq_id,
            Cartn_x: Column.ofFloatArray(x),
            Cartn_y: Column.ofFloatArray(y),
            Cartn_z: Column.ofFloatArray(z),
            id: Column.ofIntArray(seqIds.map((_, i) => i + 1)),
            label_asym_id: auth_asym_id,
            label_atom_id: auth_atom_id,
            label_comp_id: auth_comp_id,
            label_seq_id: auth_seq_id,
            label_entity_id: Column.ofConst('1', count, Column.Schema.str),
            occupancy: Column.ofConst(1, count, Column.Schema.float),
            type_symbol: Column.ofConst('C', count, Column.Schema.str),
            pdbx_PDB_model_num: Column.ofConst(1, count, Column.Schema.int),
        },
        count,
    );

    const entityBuilder = new EntityBuilder();
    for (const chain of chains) entityBuilder.getEntityId('ALA', MoleculeType.Protein, chain);
    const componentBuilder = new ComponentBuilder(auth_seq_id, auth_atom_id);
    componentBuilder.add('ALA', 0);

    const basic = createBasic({
        entity: entityBuilder.getEntityTable(),
        chem_comp: componentBuilder.getChemCompTable(),
        atom_site,
    });
    const trajectory = await createModels(
        basic,
        { kind: 'visual-test', name: 'helices', data: undefined },
        RuntimeContext.Synchronous,
    );
    return Structure.ofModel(trajectory.representative);
}

export async function structureRepr<P extends PD.Params>(
    graphics: GraphicsContext,
    provider: RepresentationProvider<Structure, P>,
    structure: Structure,
    props: Partial<PD.Values<P>>,
    colorTheme?: { name: string; params?: {} },
): Promise<Representation.Any> {
    const registries = getThemeRegistries();
    const repr = provider.factory({ graphics, ...registries }, provider.getParams);
    repr.setTheme(
        Theme.create(
            registries,
            { structure },
            {
                colorTheme: {
                    name: colorTheme?.name ?? provider.defaultColorTheme.name,
                    params: colorTheme?.params ?? provider.defaultColorTheme.props ?? {},
                },
                sizeTheme: { name: provider.defaultSizeTheme.name, params: provider.defaultSizeTheme.props ?? {} },
            },
        ),
    );
    await repr.createOrUpdate(props, structure).run();
    return repr as unknown as Representation.Any;
}
