/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit, StructureElement, StructureProperties as Props, Bond } from '../mol-model/structure';
import { Loci } from '../mol-model/loci';
import { OrderedSet } from '../mol-data/int';
import { capitalize, stripTags } from '../mol-util/string';
import { Column } from '../mol-data/db';

export type LabelGranularity = 'element' | 'conformation' | 'residue' | 'chain' | 'structure'

export const DefaultLabelOptions = {
    granularity: 'element' as LabelGranularity,
    countsOnly: false,
    hidePrefix: false,
    htmlStyling: true,
}
export type LabelOptions = typeof DefaultLabelOptions

export function lociLabel(loci: Loci, options: Partial<LabelOptions> = {}): string {
    switch (loci.kind) {
        case 'structure-loci':
            return loci.structure.models.map(m => m.entry).join(', ')
        case 'element-loci':
            return structureElementStatsLabel(StructureElement.Stats.ofLoci(loci), options)
        case 'bond-loci':
            const bond = loci.bonds[0]
            return bond ? bondLabel(bond) : ''
        case 'shape-loci':
            return loci.shape.name
        case 'group-loci':
            const g = loci.groups[0]
            return g ? loci.shape.getLabel(OrderedSet.start(g.ids), g.instance) : ''
        case 'every-loci':
            return 'Everything'
        case 'empty-loci':
            return 'Nothing'
        case 'data-loci':
            return loci.getLabel()
    }
}

function countLabel(count: number, label: string) {
    return count === 1 ? `1 ${label}` : `${count} ${label}s`
}

function otherLabel(count: number, location: StructureElement.Location, granularity: LabelGranularity, hidePrefix: boolean) {
    return `${elementLabel(location, { granularity, hidePrefix })} <small>[+ ${countLabel(count - 1, `other ${capitalize(granularity)}`)}]</small>`
}

/** Gets residue count of the model chain segments the unit is a subset of */
function getResidueCount(unit: Unit.Atomic) {
    const { elements, model } = unit
    const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy
    const elementStart = chainAtomSegments.offsets[chainAtomSegments.index[elements[0]]]
    const elementEnd = chainAtomSegments.offsets[chainAtomSegments.index[elements[elements.length - 1]] + 1]
    return residueAtomSegments.index[elementEnd] - residueAtomSegments.index[elementStart]
}

export function structureElementStatsLabel(stats: StructureElement.Stats, options: Partial<LabelOptions> = {}): string {
    const o = { ...DefaultLabelOptions, ...options }
    const label = _structureElementStatsLabel(stats, o.countsOnly, o.hidePrefix)
    return o.htmlStyling ? label : stripTags(label)
}

function _structureElementStatsLabel(stats: StructureElement.Stats, countsOnly = false, hidePrefix = false): string {
    const { structureCount, chainCount, residueCount, conformationCount, elementCount } = stats

    if (!countsOnly && elementCount === 1 && residueCount === 0 && chainCount === 0) {
        return elementLabel(stats.firstElementLoc, { hidePrefix, granularity: 'element' })
    } else if (!countsOnly && elementCount === 0 && residueCount === 1 && chainCount === 0) {
        return elementLabel(stats.firstResidueLoc, { hidePrefix, granularity: 'residue' })
    } else if (!countsOnly && elementCount === 0 && residueCount === 0 && chainCount === 1) {
        const { unit } = stats.firstChainLoc
        const granularity = (Unit.isAtomic(unit) && getResidueCount(unit) === 1) ? 'residue' : 'chain'
        return elementLabel(stats.firstChainLoc, { hidePrefix, granularity })
    } else if (!countsOnly) {
        const label: string[] = []
        if (structureCount > 0) {
            label.push(structureCount === 1 ? elementLabel(stats.firstStructureLoc, { hidePrefix, granularity: 'structure' }) : otherLabel(structureCount, stats.firstStructureLoc, 'structure', false))
        }
        if (chainCount > 0) {
            label.push(chainCount === 1 ? elementLabel(stats.firstChainLoc, { granularity: 'chain' }) : otherLabel(chainCount, stats.firstChainLoc, 'chain', false))
            hidePrefix = true;
        }
        if (residueCount > 0) {
            label.push(residueCount === 1 ? elementLabel(stats.firstResidueLoc, { granularity: 'residue', hidePrefix }) : otherLabel(residueCount, stats.firstResidueLoc, 'residue', hidePrefix))
            hidePrefix = true;
        }
        if (conformationCount > 0) {
            label.push(conformationCount === 1 ? elementLabel(stats.firstConformationLoc, { granularity: 'conformation', hidePrefix }) : otherLabel(conformationCount, stats.firstConformationLoc, 'conformation', hidePrefix))
            hidePrefix = true;
        }
        if (elementCount > 0) {
            label.push(elementCount === 1 ? elementLabel(stats.firstElementLoc, { granularity: 'element', hidePrefix }) : otherLabel(elementCount, stats.firstElementLoc, 'element', hidePrefix))
        }
        return label.join('<small> + </small>')
    } else {
        const label: string[] = []
        if (structureCount > 0) label.push(countLabel(structureCount, 'Structure'))
        if (chainCount > 0) label.push(countLabel(chainCount, 'Chain'))
        if (residueCount > 0) label.push(countLabel(residueCount, 'Residue'))
        if (conformationCount > 0) label.push(countLabel(conformationCount, 'Conformation'))
        if (elementCount > 0) label.push(countLabel(elementCount, 'Element'))
        return label.join('<small> + </small>')
    }
}

export function bondLabel(bond: Bond.Location, options: Partial<LabelOptions> = {}): string {
    const o = { ...DefaultLabelOptions, ...options }
    const locA = StructureElement.Location.create(bond.aUnit, bond.aUnit.elements[bond.aIndex])
    const locB = StructureElement.Location.create(bond.bUnit, bond.bUnit.elements[bond.bIndex])
    const labelA = _elementLabel(locA, o.granularity, o.hidePrefix)
    const labelB = _elementLabel(locB, o.granularity, o.hidePrefix)
    let offset = 0
    for (let i = 0, il = Math.min(labelA.length, labelB.length); i < il; ++i) {
        if (labelA[i] === labelB[i]) offset += 1
        else break
    }
    return `${labelA.join(' | ')} \u2014 ${labelB.slice(offset).join(' | ')}`
}

export function bundleLabel(bundle: Loci.Bundle<any>, options: Partial<LabelOptions> = {}) {
    let isSingleElements = true
    for (const l of bundle.loci) {
        if (!StructureElement.Loci.is(l) || StructureElement.Loci.size(l) !== 1) {
            isSingleElements = false
            break
        }
    }

    if (isSingleElements) {
        const o = { ...DefaultLabelOptions, ...options }
        const locations = (bundle.loci as StructureElement.Loci[]).map(l => {
            const { unit, indices } = l.elements[0]
            return StructureElement.Location.create(unit, unit.elements[OrderedSet.start(indices)])
        })
        const labels = locations.map(l => _elementLabel(l, o.granularity, o.hidePrefix))

        let offset = 0
        for (let i = 0, il = Math.min(...labels.map(l => l.length)); i < il; ++i) {
            let areIdentical = true
            for (let j = 1, jl = labels.length; j < jl; ++j) {
                if (labels[0][i] !== labels[j][i]) {
                    areIdentical = false
                    break
                }
            }
            if (areIdentical) offset += 1
            else break
        }

        if (offset > 0) {
            const offsetLabels = [labels[0].join(' | ')]
            for (let j = 1, jl = labels.length; j < jl; ++j) {
                offsetLabels.push(labels[j].slice(offset).join(' | '))
            }
            return offsetLabels.join(' \u2014 ')
        } else {
            return labels.map(l => l.join(' | ')).join('</br>')
        }
    } else {
        return bundle.loci.map(l => lociLabel(l)).join('</br>')
    }
}

export function elementLabel(location: StructureElement.Location, options: Partial<LabelOptions> = {}): string {
    const o = { ...DefaultLabelOptions, ...options }
    const label = _elementLabel(location, o.granularity, o.hidePrefix).join(' | ')
    return o.htmlStyling ? label : stripTags(label)
}

function _elementLabel(location: StructureElement.Location, granularity: LabelGranularity = 'element', hidePrefix = false): string[] {
    const label: string[] = [];
    if (!hidePrefix) {
        label.push(`<small>${location.unit.model.entry}</small>`) // entry
        if (granularity !== 'structure') {
            label.push(`<small>Model ${location.unit.model.modelNum}</small>`) // model
            label.push(`<small>Instance ${location.unit.conformation.operator.name}</small>`) // instance
        }
    }

    if (Unit.isAtomic(location.unit)) {
        label.push(..._atomicElementLabel(location as StructureElement.Location<Unit.Atomic>, granularity))
    } else if (Unit.isCoarse(location.unit)) {
        label.push(..._coarseElementLabel(location as StructureElement.Location<Unit.Spheres | Unit.Gaussians>, granularity))
    } else {
        label.push('Unknown')
    }

    return label
}

function _atomicElementLabel(location: StructureElement.Location<Unit.Atomic>, granularity: LabelGranularity): string[] {
    const label_asym_id = Props.chain.label_asym_id(location)
    const auth_asym_id = Props.chain.auth_asym_id(location)
    const has_label_seq_id = location.unit.model.atomicHierarchy.residues.label_seq_id.valueKind(location.element) === Column.ValueKind.Present;
    const label_seq_id = Props.residue.label_seq_id(location)
    const auth_seq_id = Props.residue.auth_seq_id(location)
    const ins_code = Props.residue.pdbx_PDB_ins_code(location)
    const comp_id = Props.residue.label_comp_id(location)
    const atom_id = Props.atom.label_atom_id(location)
    const alt_id = Props.atom.label_alt_id(location)
    const occupancy = Props.atom.occupancy(location);

    const microHetCompIds = Props.residue.microheterogeneityCompIds(location)
    const compId = granularity === 'residue' && microHetCompIds.length > 1 ?
        `(${microHetCompIds.join('|')})` : comp_id

    const label: string[] = []

    switch (granularity) {
        case 'element':
            label.push(`<b>${atom_id}</b>${alt_id ? `%${alt_id}` : ''}`)
        case 'conformation':
            if (granularity === 'conformation' && alt_id) {
                label.push(`<small>Conformation</small> <b>${alt_id}</b>`)
            }
        case 'residue':
            label.push(`<b>${compId}${has_label_seq_id ? ` ${label_seq_id}` : ''}</b>${label_seq_id !== auth_seq_id ? ` <small>[auth</small> <b>${auth_seq_id}</b><small>]</small>` : ''}<b>${ins_code ? ins_code : ''}</b>`)
        case 'chain':
            if (label_asym_id === auth_asym_id) {
                label.push(`<b>${label_asym_id}</b>`)
            } else {
                if (granularity === 'chain' && Unit.Traits.is(location.unit.traits, Unit.Trait.MultiChain)) {
                    label.push(`<small>[auth</small> <b>${auth_asym_id}</b><small>]</small>`)
                } else {
                    label.push(`<b>${label_asym_id}</b> <small>[auth</small> <b>${auth_asym_id}</b><small>]</small>`)
                }
            }
    }

    if (label.length > 0 && occupancy !== 1) {
        label[0] = `${label[0]} <small>[occupancy</small> <b>${Math.round(100 * occupancy) / 100}</b><small>]</small>`;
    }

    return label.reverse()
}

function _coarseElementLabel(location: StructureElement.Location<Unit.Spheres | Unit.Gaussians>, granularity: LabelGranularity): string[] {
    const asym_id = Props.coarse.asym_id(location)
    const seq_id_begin = Props.coarse.seq_id_begin(location)
    const seq_id_end = Props.coarse.seq_id_end(location)

    const label: string[] = []

    switch (granularity) {
        case 'element':
        case 'conformation':
        case 'residue':
            if (seq_id_begin === seq_id_end) {
                const entityIndex = Props.coarse.entityKey(location)
                const seq = location.unit.model.sequence.byEntityKey[entityIndex]
                const comp_id = seq.sequence.compId.value(seq_id_begin - 1) // 1-indexed
                label.push(`<b>${comp_id} ${seq_id_begin}-${seq_id_end}</b>`)
            } else {
                label.push(`<b>${seq_id_begin}-${seq_id_end}</b>`)
            }
        case 'chain':
            label.push(`<b>${asym_id}</b>`)
    }

    return label.reverse()
}