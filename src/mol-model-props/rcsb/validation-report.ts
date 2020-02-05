/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { CustomPropertyDescriptor, Structure, Unit } from '../../mol-model/structure';
// import { Database as _Database } from '../../mol-data/db'
import { CustomProperty } from '../common/custom-property';
import { CustomModelProperty } from '../common/custom-model-property';
import { Model, ElementIndex, ResidueIndex } from '../../mol-model/structure/model';
import { IntAdjacencyGraph } from '../../mol-math/graph';
import { readFromFile } from '../../mol-util/data-source';
import { CustomStructureProperty } from '../common/custom-structure-property';
import { InterUnitGraph } from '../../mol-math/graph/inter-unit-graph';
import { UnitIndex } from '../../mol-model/structure/structure/element/element';
import { IntMap, SortedArray } from '../../mol-data/int';
import { arrayMax } from '../../mol-util/array';
import { equalEps } from '../../mol-math/linear-algebra/3d/common';
import { Vec3 } from '../../mol-math/linear-algebra';

export { ValidationReport }

interface ValidationReport {
    /**
     * Real Space R (RSRZ) for residues,
     * defined for polymer residues in X-ray structures
     */
    rsrz: Map<ResidueIndex, number>
    /**
     * Real Space Correlation Coefficient (RSCC) for residues,
     * defined for each non-polymer residue in X-ray structures
     */
    rscc: Map<ResidueIndex, number>
    /**
     * Random Coil Index (RCI) for residues,
     * defined for polymer residues in NMR structures
     */
    rci: Map<ResidueIndex, number>
    /**
     * Set of geometry issues for residues
     */
    geometryIssues: Map<ResidueIndex, Set<string>>

    /**
     * Set of bond outliers
     */
    bondOutliers: {
        index: Map<ElementIndex, number>
        data: {
            atomA: ElementIndex, atomB: ElementIndex
            zScore: number, mean: number, mindiff: number,
            numobs: number, obsval: number, stdev: number
        }[]
    }
    /**
     * Set of angle outliers
     */
    angleOutliers: {
        index: Map<ElementIndex, number>
        data: {
            atomA: ElementIndex, atomB: ElementIndex, atomC: ElementIndex,
            zScore: number, mean: number, mindiff: number,
            numobs: number, obsval: number, stdev: number
        }[]
    }

    /**
     * Clashes between atoms, including id, magniture and distance
     */
    clashes: IntAdjacencyGraph<ElementIndex, {
        readonly id: ArrayLike<number>
        readonly magnitude: ArrayLike<number>
        readonly distance: ArrayLike<number>
    }>
}

namespace ValidationReport {
    export enum Tag {
        DensityFit = 'rcsb-density-fit',
        GeometryQuality = 'rcsb-geometry-quality',
        RandomCoilIndex = 'rcsb-random-coil-index',
        Clashes = 'rcsb-clashes',
    }

    export const DefaultBaseUrl = '//ftp.rcsb.org/pub/pdb/validation_reports'
    export function getEntryUrl(pdbId: string, baseUrl: string) {
        const id = pdbId.toLowerCase()
        return `${baseUrl}/${id.substr(1, 2)}/${id}/${id}_validation.xml.gz`
    }

    export function isApplicable(model?: Model): boolean {
        return (
            !!model &&
            model.sourceData.kind === 'mmCIF' &&
            (model.sourceData.data.database_2.database_id.isDefined ||
                model.entryId.length === 4)
        )
    }

    export function fromXml(xml: XMLDocument, model: Model): ValidationReport {
        return parseValidationReportXml(xml, model)
    }

    export async function fetch(ctx: CustomProperty.Context, model: Model, props: ServerSourceProps): Promise<ValidationReport> {
        const url = getEntryUrl(model.entryId, props.baseUrl)
        const xml = await ctx.fetch({ url, type: 'xml' }).runInContext(ctx.runtime)
        return fromXml(xml, model)
    }

    export async function open(ctx: CustomProperty.Context, model: Model, props: FileSourceProps): Promise<ValidationReport> {
        const xml = await readFromFile(props.input, 'xml').runInContext(ctx.runtime)
        return fromXml(xml, model)
    }

    export async function obtain(ctx: CustomProperty.Context, model: Model, props: ValidationReportProps): Promise<ValidationReport> {
        switch(props.source.name) {
            case 'file': return open(ctx, model, props.source.params)
            case 'server': return fetch(ctx, model, props.source.params)
        }
    }
}

const FileSourceParams = {
    input: PD.File({ accept: '.xml,.gz,.zip' })
}
type FileSourceProps = PD.Values<typeof FileSourceParams>

const ServerSourceParams = {
    baseUrl: PD.Text(ValidationReport.DefaultBaseUrl, { description: 'Base URL to directory tree' })
}
type ServerSourceProps = PD.Values<typeof ServerSourceParams>

export const ValidationReportParams = {
    source: PD.MappedStatic('server', {
        'file': PD.Group(FileSourceParams, { label: 'File', isFlat: true }),
        'server': PD.Group(ServerSourceParams, { label: 'Server', isFlat: true }),
    }, { options: [['file', 'File'], ['server', 'Server']] })
}
export type ValidationReportParams = typeof ValidationReportParams
export type ValidationReportProps = PD.Values<ValidationReportParams>

export const ValidationReportProvider: CustomModelProperty.Provider<ValidationReportParams, ValidationReport> = CustomModelProperty.createProvider({
    label: 'RCSB Validation Report',
    descriptor: CustomPropertyDescriptor({
        name: 'rcsb_validation_report',
        // TODO `cifExport` and `symbol`
    }),
    type: 'dynamic',
    defaultParams: ValidationReportParams,
    getParams: (data: Model) => ValidationReportParams,
    isApplicable: (data: Model) => ValidationReport.isApplicable(data),
    obtain: async (ctx: CustomProperty.Context, data: Model, props: Partial<ValidationReportProps>) => {
        const p = { ...PD.getDefaultValues(ValidationReportParams), ...props }
        return await ValidationReport.obtain(ctx, data, p)
    }
})

//

type IntraUnitClashesProps = {
    readonly id: ArrayLike<number>
    readonly magnitude: ArrayLike<number>
    readonly distance: ArrayLike<number>
}
type InterUnitClashesProps = {
    readonly id: number
    readonly magnitude: number
    readonly distance: number
}

export type IntraUnitClashes = IntAdjacencyGraph<UnitIndex, IntraUnitClashesProps>
export type InterUnitClashes = InterUnitGraph<Unit.Atomic, UnitIndex, InterUnitClashesProps>

export interface Clashes {
    readonly interUnit: InterUnitClashes
    readonly intraUnit: IntMap<IntraUnitClashes>
}

function createInterUnitClashes(structure: Structure, clashes: ValidationReport['clashes']) {
    const builder = new InterUnitGraph.Builder<Unit.Atomic, UnitIndex, InterUnitClashesProps>()
    const { a, b, edgeProps: { id, magnitude, distance } } = clashes

    const pA = Vec3()
    const pB = Vec3()

    Structure.eachUnitPair(structure, (unitA: Unit, unitB: Unit) => {
        const elementsA = unitA.elements
        const elementsB = unitB.elements

        builder.startUnitPair(unitA as Unit.Atomic, unitB as Unit.Atomic)

        for (let i = 0, il = clashes.edgeCount * 2; i < il; ++i) {
            // TODO create lookup
            let indexA = SortedArray.indexOf(elementsA, a[i])
            let indexB = SortedArray.indexOf(elementsB, b[i])

            if (indexA !== -1 && indexB !== -1) {
                unitA.conformation.position(a[i], pA)
                unitB.conformation.position(b[i], pB)

                // check actual distance to avoid clashes between unrelated chain instances
                if (equalEps(distance[i], Vec3.distance(pA, pB), 0.1)) {
                    builder.add(indexA as UnitIndex, indexB as UnitIndex, {
                        id: id[i],
                        magnitude: magnitude[i],
                        distance: distance[i]
                    })
                }
            }
        }

        builder.finishUnitPair()
    }, {
        maxRadius: arrayMax(clashes.edgeProps.distance),
        validUnit: (unit: Unit) => Unit.isAtomic(unit),
        validUnitPair: (unitA: Unit, unitB: Unit) => unitA.model === unitB.model
    })

    return new InterUnitGraph(builder.getMap())
}

function createIntraUnitClashes(unit: Unit.Atomic, clashes: ValidationReport['clashes']): IntraUnitClashes {
    const aIndices: UnitIndex[] = []
    const bIndices: UnitIndex[] = []
    const ids: number[] = []
    const magnitudes: number[] = []
    const distances: number[] = []

    const { elements } = unit
    const { a, b, edgeCount, edgeProps } = clashes

    for (let i = 0, il = edgeCount * 2; i < il; ++i) {
        // TODO create lookup
        let indexA = SortedArray.indexOf(elements, a[i])
        let indexB = SortedArray.indexOf(elements, b[i])

        if (indexA !== -1 && indexB !== -1) {
            aIndices.push(indexA as UnitIndex)
            bIndices.push(indexB as UnitIndex)
            ids.push(edgeProps.id[i])
            magnitudes.push(edgeProps.magnitude[i])
            distances.push(edgeProps.distance[i])
        }
    }

    const builder = new IntAdjacencyGraph.EdgeBuilder(elements.length, aIndices, bIndices)
    const id = new Int32Array(builder.slotCount)
    const magnitude = new Float32Array(builder.slotCount)
    const distance = new Float32Array(builder.slotCount)
    for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
        builder.addNextEdge()
        builder.assignProperty(id, ids[i])
        builder.assignProperty(magnitude, magnitudes[i])
        builder.assignProperty(distance, distances[i])
    }
    return builder.createGraph({ id, magnitude, distance })
}

function createClashes(structure: Structure, clashes: ValidationReport['clashes']): Clashes {

    const intraUnit = IntMap.Mutable<IntraUnitClashes>()

    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
        const group = structure.unitSymmetryGroups[i]
        if (!Unit.isAtomic(group.units[0])) continue

        const intraClashes = createIntraUnitClashes(group.units[0], clashes)
        for (let j = 0, jl = group.units.length; j < jl; ++j) {
            intraUnit.set(group.units[j].id, intraClashes)
        }
    }

    return {
        interUnit: createInterUnitClashes(structure, clashes),
        intraUnit
    }
}

export const ClashesProvider: CustomStructureProperty.Provider<{}, Clashes> = CustomStructureProperty.createProvider({
    label: 'Clashes',
    descriptor: CustomPropertyDescriptor({
        name: 'rcsb_clashes',
        // TODO `cifExport` and `symbol`
    }),
    type: 'local',
    defaultParams: {},
    getParams: (data: Structure) => ({}),
    isApplicable: (data: Structure) => true,
    obtain: async (ctx: CustomProperty.Context, data: Structure) => {
        await ValidationReportProvider.attach(ctx, data.models[0])
        const validationReport = ValidationReportProvider.get(data.models[0]).value!
        return createClashes(data, validationReport.clashes)
    }
})

//

function getItem(a: NamedNodeMap, name: string) {
    const item = a.getNamedItem(name)
    return item !== null ? item.value : ''
}

function hasAttr(a: NamedNodeMap, name: string, value: string) {
    const item = a.getNamedItem(name)
    return item !== null && item.value === value
}

function getMogInfo(a: NamedNodeMap) {
    return {
        zScore: parseFloat(getItem(a, 'Zscore')),
        mean: parseFloat(getItem(a, 'mean')),
        mindiff: parseFloat(getItem(a, 'mindiff')),
        numobs: parseInt(getItem(a, 'numobs')),
        obsval: parseFloat(getItem(a, 'obsval')),
        stdev: parseFloat(getItem(a, 'stdev')),
    }
}

function ClashesBuilder(elementsCount: number) {
    const aIndices: ElementIndex[] = []
    const bIndices: ElementIndex[] = []
    const ids: number[] = []
    const magnitudes: number[] = []
    const distances: number[] = []

    const seen = new Map<number, ElementIndex>()

    return {
        add(element: ElementIndex, id: number, magnitude: number, distance: number) {
            const other = seen.get(id)
            if (other !== undefined) {
                aIndices[aIndices.length] = element
                bIndices[bIndices.length] = other
                ids[ids.length] = id
                magnitudes[magnitudes.length] = magnitude
                distances[distances.length] = distance
            } else {
                seen.set(id, element)
            }
        },
        get() {
            const builder = new IntAdjacencyGraph.EdgeBuilder(elementsCount, aIndices, bIndices)
            const id = new Int32Array(builder.slotCount)
            const magnitude = new Float32Array(builder.slotCount)
            const distance = new Float32Array(builder.slotCount)
            for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
                builder.addNextEdge()
                builder.assignProperty(id, ids[i])
                builder.assignProperty(magnitude, magnitudes[i])
                builder.assignProperty(distance, distances[i])
            }
            return builder.createGraph({ id, magnitude, distance })
        }
    }
}

function parseValidationReportXml(xml: XMLDocument, model: Model): ValidationReport {
    const rsrz = new Map<ResidueIndex, number>()
    const rscc = new Map<ResidueIndex, number>()
    const rci = new Map<ResidueIndex, number>()
    const geometryIssues = new Map<ResidueIndex, Set<string>>()

    const bondOutliers = {
        index: new Map<ElementIndex, number>(),
        data: [] as ValidationReport['bondOutliers']['data']
    }
    const angleOutliers = {
        index: new Map<ElementIndex, number>(),
        data: [] as ValidationReport['angleOutliers']['data']
    }

    const clashesBuilder = ClashesBuilder(model.atomicHierarchy.atoms._rowCount)

    const { index } = model.atomicHierarchy

    const entries = xml.getElementsByTagName('Entry')
    if (entries.length === 1) {
        const chemicalShiftLists = entries[0].getElementsByTagName('chemical_shift_list')
        if (chemicalShiftLists.length === 1) {
            const randomCoilIndices = chemicalShiftLists[0].getElementsByTagName('random_coil_index')
            for (let j = 0, jl = randomCoilIndices.length; j < jl; ++j) {
                const { attributes } = randomCoilIndices[j]
                const value = parseFloat(getItem(attributes, 'value'))
                const auth_asym_id = getItem(attributes, 'chain')
                const auth_comp_id = getItem(attributes, 'rescode')
                const auth_seq_id = parseInt(getItem(attributes, 'resnum'))
                const rI = index.findResidueAuth({ auth_asym_id, auth_comp_id, auth_seq_id })
                if (rI !== -1) rci.set(rI, value)
            }
        }
    }

    const groups = xml.getElementsByTagName('ModelledSubgroup')
    for (let i = 0, il = groups.length; i < il; ++i) {
        const g = groups[ i ]
        const ga = g.attributes

        const pdbx_PDB_model_num = parseInt(getItem(ga, 'model'))
        if (model.modelNum !== pdbx_PDB_model_num) continue

        const auth_asym_id = getItem(ga, 'chain')
        const auth_comp_id = getItem(ga, 'resname')
        const auth_seq_id = parseInt(getItem(ga, 'resnum'))
        const pdbx_PDB_ins_code = getItem(ga, 'icode').trim() || undefined
        const label_alt_id = getItem(ga, 'altcode').trim() || undefined

        const rI = index.findResidueAuth({ auth_asym_id, auth_comp_id, auth_seq_id, pdbx_PDB_ins_code })

        // continue if no residue index is found
        if (rI === -1) continue

        if (ga.getNamedItem('rsrz') !== null) rsrz.set(rI, parseFloat(getItem(ga, 'rsrz')))
        if (ga.getNamedItem('rscc') !== null) rscc.set(rI, parseFloat(getItem(ga, 'rscc')))

        const isPolymer = getItem(ga, 'seq') !== '.'
        const issues = new Set<string>()

        if (isPolymer) {
            const angleOutliers = g.getElementsByTagName('angle-outlier')
            if (angleOutliers.length) issues.add('angle-outlier')

            const bondOutliers = g.getElementsByTagName('bond-outlier')
            if (bondOutliers.length) issues.add('bond-outlier')

            const planeOutliers = g.getElementsByTagName('plane-outlier')
            if (planeOutliers.length) issues.add('plane-outlier')

            if (hasAttr(ga, 'rota', 'OUTLIER')) issues.add('rotamer-outlier')
            if (hasAttr(ga, 'rama', 'OUTLIER')) issues.add('ramachandran-outlier')
            if (hasAttr(ga, 'RNApucker', 'outlier')) issues.add('RNApucker-outlier')
        } else {
            const mogBondOutliers = g.getElementsByTagName('mog-bond-outlier')
            if (mogBondOutliers.length) issues.add('mogul-bond-outlier')

            const mogAngleOutliers = g.getElementsByTagName('mog-angle-outlier')
            if (mogAngleOutliers.length) issues.add('mogul-angle-outlier')

            for (let j = 0, jl = mogBondOutliers.length; j < jl; ++j) {
                const mbo = mogBondOutliers[ j ].attributes
                const atoms = getItem(mbo, 'atoms').split(',')
                const idx = bondOutliers.data.length
                const atomA = index.findAtomOnResidue(rI, atoms[0])
                const atomB = index.findAtomOnResidue(rI, atoms[1])
                bondOutliers.index.set(atomA, idx)
                bondOutliers.index.set(atomB, idx)
                bondOutliers.data.push({ atomA, atomB, ...getMogInfo(mbo) })
            }

            for (let j = 0, jl = mogAngleOutliers.length; j < jl; ++j) {
                const mao = mogAngleOutliers[ j ].attributes
                const atoms = getItem(mao, 'atoms').split(',')
                const idx = angleOutliers.data.length
                const atomA = index.findAtomOnResidue(rI, atoms[0])
                const atomB = index.findAtomOnResidue(rI, atoms[1])
                const atomC = index.findAtomOnResidue(rI, atoms[1])
                angleOutliers.index.set(atomA, idx)
                angleOutliers.index.set(atomB, idx)
                angleOutliers.index.set(atomC, idx)
                angleOutliers.data.push({ atomA, atomB, atomC, ...getMogInfo(mao) })
            }
        }

        const clashes = g.getElementsByTagName('clash')
        if (clashes.length) issues.add('clash')

        for (let j = 0, jl = clashes.length; j < jl; ++j) {
            const ca = clashes[ j ].attributes
            const id = parseInt(getItem(ca, 'cid'))
            const magnitude = parseFloat(getItem(ca, 'clashmag'))
            const distance = parseFloat(getItem(ca, 'dist'))
            const label_atom_id = getItem(ca, 'atom')
            const element = index.findAtomOnResidue(rI, label_atom_id, label_alt_id)
            if (element !== -1) {
                clashesBuilder.add(element, id, magnitude, distance)
            }
        }

        geometryIssues.set(rI, issues)
    }

    const clashes = clashesBuilder.get()

    const validationReport = {
        rsrz, rscc, rci, geometryIssues,
        bondOutliers, angleOutliers,
        clashes
    }
    console.log(validationReport)

    return validationReport
}