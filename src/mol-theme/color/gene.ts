/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureProperties, StructureElement, Link } from 'mol-model/structure';
import { ColorScale, Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from 'mol-theme/theme';
import { ColorListOptions, ColorListName } from 'mol-util/color/scale';
import { NumberArray } from 'mol-util/type-helpers';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives ranges of a polymer chain a color based on the gene (or linker/terminal extension) it originates from.'

export const GeneColorThemeParams = {
    list: PD.ColorScale<ColorListName>('RedYellowBlue', ColorListOptions),
}
export type GeneColorThemeParams = typeof GeneColorThemeParams
export function getGeneColorThemeParams(ctx: ThemeDataContext) {
    return GeneColorThemeParams // TODO return copy
}

function modelEntityKey(modelIndex: number, entityId: string) {
    return `${modelIndex}|${entityId}`
}

function addGene(geneSerialMap: Map<string, number>, geneNames: string[], beg: number, end: number, seqToSrcGen: NumberArray) {
    const gene = geneNames.map(s => s.toUpperCase()).sort().join(',')
    let geneIndex = 0 // serial no starting from 1
    if (gene === '') {
        geneIndex = geneSerialMap.size + 1
        geneSerialMap.set(`UNKNOWN${geneIndex}`, geneIndex)
    } else if (geneSerialMap.has(gene)) {
        geneIndex = geneSerialMap.get(gene)!
    } else {
        geneIndex = geneSerialMap.size + 1
        geneSerialMap.set(gene, geneIndex)
    }
    for (let i = beg, il = end; i <= il; ++i) {
        seqToSrcGen[i - 1] = geneIndex
    }
}

export function GeneColorTheme(ctx: ThemeDataContext, props: PD.Values<GeneColorThemeParams>): ColorTheme<GeneColorThemeParams> {
    let color: LocationColor
    const scale = ColorScale.create({ listOrName: props.list, minLabel: 'Start', maxLabel: 'End' })
    const { structure } = ctx

    if (structure) {
        const l = StructureElement.create()
        const { models } = structure
        const seqToSrcGenByModelEntity = new Map<string, NumberArray>()
        const geneSerialMap = new Map<string, number>() // serial no starting from 1
        for (let i = 0, il = models.length; i <il; ++i) {
            const m = models[i]
            if (m.sourceData.kind !== 'mmCIF') continue
            const { entity_src_gen } = m.sourceData.data
            
            const { entity_id, pdbx_beg_seq_num, pdbx_end_seq_num, pdbx_gene_src_gene } = entity_src_gen
            for (let j = 0, jl = entity_src_gen._rowCount; j < jl; ++j) {
                const entityId = entity_id.value(j)
                const k = modelEntityKey(i, entityId)
                if (!seqToSrcGenByModelEntity.has(k)) {
                    const entityIndex = m.entities.getEntityIndex(entityId)
                    const seq = m.sequence.sequences[entityIndex].sequence
                    const seqLength = seq.sequence.length
                    const seqToGene = new Int16Array(seqLength)
                    addGene(geneSerialMap, pdbx_gene_src_gene.value(j), pdbx_beg_seq_num.value(j), pdbx_end_seq_num.value(j), seqToGene)
                    seqToSrcGenByModelEntity.set(k, seqToGene)
                } else {
                    const seqToGene = seqToSrcGenByModelEntity.get(k)!
                    addGene(geneSerialMap, pdbx_gene_src_gene.value(j), pdbx_beg_seq_num.value(j), pdbx_end_seq_num.value(j), seqToGene)
                    seqToSrcGenByModelEntity.set(k, seqToGene)
                }
            }
        }
        scale.setDomain(1, geneSerialMap.size)
        const scaleColor = scale.color

        const getGeneColor = (location: StructureElement) => {
            const modelIndex = structure.models.indexOf(location.unit.model)
            const entityId = StructureProperties.entity.id(location)
            const k = modelEntityKey(modelIndex, entityId)
            const seqToGene = seqToSrcGenByModelEntity.get(k)
            if (seqToGene) {
                // minus 1 to convert seqId to array index
                return scaleColor(seqToGene[StructureProperties.residue.label_seq_id(location) - 1])
            } else {
                return DefaultColor
            }
        }

        color = (location: Location): Color => {
            if (StructureElement.isLocation(location)) {
                return getGeneColor(location)
            } else if (Link.isLocation(location)) {
                l.unit = location.aUnit
                l.element = location.aUnit.elements[location.aIndex]
                return getGeneColor(l)
            }
            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        factory: GeneColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: scale ? scale.legend : undefined
    }
}

export const GeneColorThemeProvider: ColorTheme.Provider<GeneColorThemeParams> = {
    label: 'Gene',
    factory: GeneColorTheme,
    getParams: getGeneColorThemeParams,
    defaultValues: PD.getDefaultValues(GeneColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
}