/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureProperties, StructureElement, Bond, Model } from '../../mol-model/structure';
import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { Table, Column } from '../../mol-data/db';
import { mmCIF_Schema } from '../../mol-io/reader/cif/schema/mmcif';
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { TableLegend, ScaleLegend } from '../../mol-util/legend';
import { isInteger } from '../../mol-util/number';
import { ColorLists, getColorListFromName } from '../../mol-util/color/lists';
import { MmcifFormat } from '../../mol-model-formats/structure/mmcif';

const DefaultList = 'dark-2';
const DefaultColor = Color(0xFAFAFA);
const Description = 'Gives ranges of a polymer chain a color based on the entity source it originates from (e.g. gene, plasmid, organism).';

export const EntitySourceColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type EntitySourceColorThemeParams = typeof EntitySourceColorThemeParams
export function getEntitySourceColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(EntitySourceColorThemeParams);
    if (ctx.structure) {
        if (getMaps(ctx.structure.root.models).srcKeySerialMap.size > ColorLists[DefaultList].list.length) {
            params.palette.defaultValue.name = 'colors';
            params.palette.defaultValue.params = {
                ...params.palette.defaultValue.params,
                list: { kind: 'interpolate', colors: getColorListFromName(DefaultList).list }
            };
        }
    }
    return params;
}

function modelEntityKey(modelIndex: number, entityId: string) {
    return `${modelIndex}|${entityId}`;
}

type EntitySrc = Table<{
    entity_id: mmCIF_Schema['entity_src_gen']['entity_id'],
    pdbx_src_id: mmCIF_Schema['entity_src_gen']['pdbx_src_id'],
    pdbx_beg_seq_num: mmCIF_Schema['entity_src_gen']['pdbx_beg_seq_num'],
    pdbx_end_seq_num: mmCIF_Schema['entity_src_gen']['pdbx_end_seq_num'],
}>
type ScientificName = Column<mmCIF_Schema['entity_src_gen']['pdbx_gene_src_scientific_name']['T']>
type GeneSrcGene = Column<mmCIF_Schema['entity_src_gen']['pdbx_gene_src_gene']['T']>
type PlasmidName = Column<mmCIF_Schema['entity_src_gen']['plasmid_name']['T']>

function srcKey(modelIndex: number, entityId: string, organism: string, srcId: number, plasmid: string, gene: string) {
    return `${modelIndex}|${entityId}|${organism}|${gene ? gene : (plasmid ? plasmid : srcId)}`;
}

function addSrc(seqToSrcByModelEntity: Map<string, Int16Array>, srcKeySerialMap: Map<string, number>, modelIndex: number, model: Model, entity_src: EntitySrc, scientific_name: ScientificName, plasmid_name?: PlasmidName, gene_src_gene?: GeneSrcGene) {
    const { entity_id, pdbx_src_id, pdbx_beg_seq_num, pdbx_end_seq_num } = entity_src;
    for (let j = 0, jl = entity_src._rowCount; j < jl; ++j) {
        const entityId = entity_id.value(j);
        const mK = modelEntityKey(modelIndex, entityId);
        let seqToSrc: Int16Array;
        if (!seqToSrcByModelEntity.has(mK)) {
            const entityIndex = model.entities.getEntityIndex(entityId);
            const seq = model.sequence.sequences[entityIndex].sequence;
            seqToSrc = new Int16Array(seq.length);
            seqToSrcByModelEntity.set(mK, seqToSrc);
        } else {
            seqToSrc = seqToSrcByModelEntity.get(mK)!;
        }

        const plasmid = plasmid_name ? plasmid_name.value(j) : '';
        const gene = gene_src_gene ? gene_src_gene.value(j)[0] : '';
        const sK = srcKey(modelIndex, entityId, scientific_name.value(j), pdbx_src_id.value(j), plasmid, gene);

        // may not be given (= 0) indicating src is for the whole seq
        const beg = pdbx_beg_seq_num.valueKind(j) === Column.ValueKind.Present ? pdbx_beg_seq_num.value(j) : 1;
        const end = pdbx_end_seq_num.valueKind(j) === Column.ValueKind.Present ? pdbx_end_seq_num.value(j) : seqToSrc.length;

        let srcIndex: number; // serial no starting from 1
        if (srcKeySerialMap.has(sK)) {
            srcIndex = srcKeySerialMap.get(sK)!;
        } else {
            srcIndex = srcKeySerialMap.size + 1;
            srcKeySerialMap.set(sK, srcIndex);
        }
        // set src index
        for (let i = beg, il = end; i <= il; ++i) {
            seqToSrc[i - 1] = srcIndex;
        }
    }
}

function getMaps(models: ReadonlyArray<Model>) {
    const seqToSrcByModelEntity = new Map<string, Int16Array>();
    const srcKeySerialMap = new Map<string, number>(); // serial no starting from 1

    for (let i = 0, il = models.length; i < il; ++i) {
        const m = models[i];
        if (!MmcifFormat.is(m.sourceData)) continue;
        const { entity_src_gen, entity_src_nat, pdbx_entity_src_syn } = m.sourceData.data.db;
        addSrc(seqToSrcByModelEntity, srcKeySerialMap, i, m, entity_src_gen, entity_src_gen.pdbx_gene_src_scientific_name, entity_src_gen.plasmid_name, entity_src_gen.pdbx_gene_src_gene);
        addSrc(seqToSrcByModelEntity, srcKeySerialMap, i, m, entity_src_nat, entity_src_nat.pdbx_organism_scientific, entity_src_nat.pdbx_plasmid_name);
        addSrc(seqToSrcByModelEntity, srcKeySerialMap, i, m, pdbx_entity_src_syn, pdbx_entity_src_syn.organism_scientific);
    }

    return { seqToSrcByModelEntity, srcKeySerialMap };
}

function getLabelTable(srcKeySerialMap: Map<string, number>) {
    let unnamedCount = 0;
    return Array.from(srcKeySerialMap.keys()).map(v => {
        const vs = v.split('|');
        const organism = vs[2];
        const name = isInteger(vs[3]) ? `Unnamed ${++unnamedCount}` : vs[3];
        return `${name}${organism ? ` (${organism})` : ''}`;
    });
}

export function EntitySourceColorTheme(ctx: ThemeDataContext, props: PD.Values<EntitySourceColorThemeParams>): ColorTheme<EntitySourceColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    if (ctx.structure) {
        const l = StructureElement.Location.create(ctx.structure);
        const { models } = ctx.structure.root;
        const { seqToSrcByModelEntity, srcKeySerialMap } = getMaps(models);

        const labelTable = getLabelTable(srcKeySerialMap);
        props.palette.params.valueLabel = (i: number) => labelTable[i];

        const palette = getPalette(srcKeySerialMap.size, props);
        legend = palette.legend;

        const getSrcColor = (location: StructureElement.Location) => {
            const modelIndex = models.indexOf(location.unit.model);
            const entityId = StructureProperties.entity.id(location);
            const mK = modelEntityKey(modelIndex, entityId);
            const seqToSrc = seqToSrcByModelEntity.get(mK);
            if (seqToSrc) {
                // minus 1 to convert seqId to array index
                return palette.color(seqToSrc[StructureProperties.residue.label_seq_id(location) - 1] - 1);
            } else {
                return DefaultColor;
            }
        };

        color = (location: Location): Color => {
            if (StructureElement.Location.is(location)) {
                return getSrcColor(location);
            } else if (Bond.isLocation(location)) {
                l.unit = location.aUnit;
                l.element = location.aUnit.elements[location.aIndex];
                return getSrcColor(l);
            }
            return DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: EntitySourceColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend
    };
}

export const EntitySourceColorThemeProvider: ColorTheme.Provider<EntitySourceColorThemeParams, 'entity-source'> = {
    name: 'entity-source',
    label: 'Entity Source',
    category: ColorTheme.Category.Chain,
    factory: EntitySourceColorTheme,
    getParams: getEntitySourceColorThemeParams,
    defaultValues: PD.getDefaultValues(EntitySourceColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};