/**
 * Copyright (c) 2017-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Segmentation, SortedArray } from '../../../../mol-data/int';
import { Structure, Unit } from '../../structure';
import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';
import { UniqueStructuresBuilder } from '../utils/builders';
import { StructureUniqueSubsetBuilder } from '../../structure/util/unique-subset-builder';
import { QueryContext, QueryFn } from '../context';
import { structureIntersect, structureSubtract, structureUnion } from '../utils/structure-set';
import { UniqueArray } from '../../../../mol-data/generic';
import { StructureSubsetBuilder } from '../../structure/util/subset-builder';
import { StructureElement } from '../../structure/element';
import { MmcifFormat } from '../../../../mol-model-formats/structure/mmcif';
import { ResidueSet, ResidueSetEntry } from '../../model/properties/utils/residue-set';
import { StructureProperties } from '../../structure/properties';
import { arraySetAdd } from '../../../../mol-util/array';

function getWholeResidues(ctx: QueryContext, source: Structure, structure: Structure) {
    const builder = source.subsetBuilder(true);
    for (const unit of structure.units) {
        if (unit.kind !== Unit.Kind.Atomic) {
            // just copy non-atomic units.
            builder.setUnit(unit.id, unit.elements);
            continue;
        }

        const { residueAtomSegments } = unit.model.atomicHierarchy;
        const sourceElements = source.unitMap.get(unit.id).elements;

        const elements = unit.elements;
        builder.beginUnit(unit.id);
        const residuesIt = Segmentation.transientSegments(residueAtomSegments, elements);
        while (residuesIt.hasNext) {
            const rI = residuesIt.move().index;
            for (let j = residueAtomSegments.offsets[rI], _j = residueAtomSegments.offsets[rI + 1]; j < _j; j++) {
                if (SortedArray.has(sourceElements, j)) builder.addElement(j);
            }
        }
        builder.commitUnit();

        ctx.throwIfTimedOut();
    }
    return builder.getStructure();
}

export function wholeResidues(query: StructureQuery): StructureQuery {
    return function query_wholeResidues(ctx) {
        const inner = query(ctx);
        if (StructureSelection.isSingleton(inner)) {
            return StructureSelection.Singletons(ctx.inputStructure, getWholeResidues(ctx, ctx.inputStructure, inner.structure));
        } else {
            const builder = new UniqueStructuresBuilder(ctx.inputStructure);
            for (const s of inner.structures) {
                builder.add(getWholeResidues(ctx, ctx.inputStructure, s));
            }
            return builder.getSelection();
        }
    };
}

export interface IncludeSurroundingsParams {
    radius: number,
    elementRadius?: QueryFn<number>,
    wholeResidues?: boolean
}

function getIncludeSurroundings(ctx: QueryContext, source: Structure, structure: Structure, params: IncludeSurroundingsParams) {
    const builder = new StructureUniqueSubsetBuilder(source);
    const lookup = source.lookup3d;
    const r = params.radius;

    for (const unit of structure.units) {
        const c = unit.conformation;
        const elements = unit.elements;
        for (let i = 0, _i = elements.length; i < _i; i++) {
            const e = elements[i];
            lookup.findIntoBuilder(c.x(e), c.y(e), c.z(e), r, builder);
        }

        ctx.throwIfTimedOut();
    }
    return !!params.wholeResidues ? getWholeResidues(ctx, source, builder.getStructure()) : builder.getStructure();
}

interface IncludeSurroundingsParamsWithRadius extends IncludeSurroundingsParams {
    elementRadius: QueryFn<number>,
    elementRadiusClosure: StructureElement.Property<number>,
    sourceMaxRadius: number
}

function getIncludeSurroundingsWithRadius(ctx: QueryContext, source: Structure, structure: Structure, params: IncludeSurroundingsParamsWithRadius) {
    const builder = new StructureUniqueSubsetBuilder(source);
    const lookup = source.lookup3d;
    const { elementRadius, elementRadiusClosure, sourceMaxRadius, radius } = params;

    ctx.pushCurrentElement();
    ctx.element.structure = structure;
    for (const unit of structure.units) {
        ctx.element.unit = unit;
        const c = unit.conformation;
        const elements = unit.elements;

        for (let i = 0, _i = elements.length; i < _i; i++) {
            const e = elements[i];
            ctx.element.element = e;
            const eRadius = elementRadius(ctx);
            lookup.findIntoBuilderWithRadius(c.x(e), c.y(e), c.z(e), eRadius, sourceMaxRadius, radius, elementRadiusClosure, builder);
        }

        ctx.throwIfTimedOut();
    }

    ctx.popCurrentElement();
    return !!params.wholeResidues ? getWholeResidues(ctx, source, builder.getStructure()) : builder.getStructure();
}

function createElementRadiusFn(ctx: QueryContext, eRadius: QueryFn<number>): StructureElement.Property<number> {
    return e => {
        ctx.element.structure = e.structure;
        ctx.element.unit = e.unit;
        ctx.element.element = e.element;
        return eRadius(ctx);
    };
}

function findStructureRadius(ctx: QueryContext, eRadius: QueryFn<number>) {
    let r = 0;
    ctx.element.structure = ctx.inputStructure;
    for (const unit of ctx.inputStructure.units) {
        ctx.element.unit = unit;
        const elements = unit.elements;

        for (let i = 0, _i = elements.length; i < _i; i++) {
            const e = elements[i];
            ctx.element.element = e;
            const eR = eRadius(ctx);
            if (eR > r) r = eR;
        }

    }
    ctx.throwIfTimedOut();
    return r;
}

export function includeSurroundings(query: StructureQuery, params: IncludeSurroundingsParams): StructureQuery {
    return function query_includeSurroundings(ctx) {
        const inner = query(ctx);

        if (params.elementRadius) {
            const prms: IncludeSurroundingsParamsWithRadius = {
                ...params,
                elementRadius: params.elementRadius,
                elementRadiusClosure: createElementRadiusFn(ctx, params.elementRadius),
                sourceMaxRadius: findStructureRadius(ctx, params.elementRadius)
            };

            if (StructureSelection.isSingleton(inner)) {
                const surr = getIncludeSurroundingsWithRadius(ctx, ctx.inputStructure, inner.structure, prms);
                const ret = StructureSelection.Singletons(ctx.inputStructure, surr);
                return ret;
            } else {
                const builder = new UniqueStructuresBuilder(ctx.inputStructure);
                for (const s of inner.structures) {
                    builder.add(getIncludeSurroundingsWithRadius(ctx, ctx.inputStructure, s, prms));
                }
                return builder.getSelection();
            }
        }

        if (StructureSelection.isSingleton(inner)) {
            const surr = getIncludeSurroundings(ctx, ctx.inputStructure, inner.structure, params);
            const ret = StructureSelection.Singletons(ctx.inputStructure, surr);
            return ret;
        } else {
            const builder = new UniqueStructuresBuilder(ctx.inputStructure);
            for (const s of inner.structures) {
                builder.add(getIncludeSurroundings(ctx, ctx.inputStructure, s, params));
            }
            return builder.getSelection();
        }
    };
}

export function querySelection(selection: StructureQuery, query: StructureQuery): StructureQuery {
    return function query_querySelection(ctx) {
        const targetSel = selection(ctx);
        if (StructureSelection.structureCount(targetSel) === 0) return targetSel;

        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);
        const add = (s: Structure) => ret.add(s);

        StructureSelection.forEach(targetSel, (s, sI) => {
            ctx.pushInputStructure(s);
            StructureSelection.forEach(query(ctx), add);
            ctx.popInputStructure();
            if (sI % 10 === 0) ctx.throwIfTimedOut();
        });
        return ret.getSelection();
    };
}

export function intersectBy(query: StructureQuery, by: StructureQuery): StructureQuery {
    return function query_intersectBy(ctx) {
        const selection = query(ctx);
        if (StructureSelection.structureCount(selection) === 0) return selection;

        const bySel = by(ctx);
        if (StructureSelection.structureCount(bySel) === 0) return StructureSelection.Empty(ctx.inputStructure);
        const unionBy = StructureSelection.unionStructure(bySel);

        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);
        StructureSelection.forEach(selection, (s, sI) => {
            const ii = structureIntersect(unionBy, s);
            if (ii.elementCount !== 0) ret.add(ii);
            if (sI % 50 === 0) ctx.throwIfTimedOut();
        });

        return ret.getSelection();
    };
}

export function exceptBy(query: StructureQuery, by: StructureQuery): StructureQuery {
    return function query_exceptBy(ctx) {
        const selection = query(ctx);
        if (StructureSelection.structureCount(selection) === 0) return selection;

        const bySel = by(ctx);
        if (StructureSelection.structureCount(bySel) === 0) return selection;
        const subtractBy = StructureSelection.unionStructure(bySel);

        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);
        StructureSelection.forEach(selection, (s, sI) => {
            const diff = structureSubtract(s, subtractBy);
            if (diff.elementCount !== 0) ret.add(diff);
            if (sI % 50 === 0) ctx.throwIfTimedOut();
        });

        return ret.getSelection();
    };
}

export function union(query: StructureQuery): StructureQuery {
    return function query_union(ctx) {
        const ret = StructureSelection.LinearBuilder(ctx.inputStructure);
        ret.add(StructureSelection.unionStructure(query(ctx)));
        return ret.getSelection();
    };
}

export function expandProperty(query: StructureQuery, property: QueryFn): StructureQuery {
    return function query_expandProperty(ctx) {
        const src = query(ctx);
        const propertyToStructureIndexMap = new Map<any, UniqueArray<number>>();

        const builders: StructureSubsetBuilder[] = [];
        ctx.pushCurrentElement();
        StructureSelection.forEach(src, (s, sI) => {
            ctx.element.structure = s;
            for (const unit of s.units) {
                ctx.element.unit = unit;
                const elements = unit.elements;
                for (let i = 0, _i = elements.length; i < _i; i++) {
                    ctx.element.element = elements[i];
                    const p = property(ctx);
                    let arr: UniqueArray<number>;
                    if (propertyToStructureIndexMap.has(p)) arr = propertyToStructureIndexMap.get(p)!;
                    else {
                        arr = UniqueArray.create<number>();
                        propertyToStructureIndexMap.set(p, arr);
                    }
                    UniqueArray.add(arr, sI, sI);
                }
            }
            builders[sI] = ctx.inputStructure.subsetBuilder(true);

            if (sI % 10 === 0) ctx.throwIfTimedOut();
        });

        ctx.element.structure = ctx.inputStructure;
        for (const unit of ctx.inputStructure.units) {
            ctx.element.unit = unit;
            const elements = unit.elements;
            for (let i = 0, _i = elements.length; i < _i; i++) {
                ctx.element.element = elements[i];
                const p = property(ctx);
                if (!propertyToStructureIndexMap.has(p)) continue;
                const indices = propertyToStructureIndexMap.get(p)!.array;

                for (let _sI = 0, __sI = indices.length; _sI < __sI; _sI++) {
                    builders[indices[_sI]].addToUnit(unit.id, elements[i]);
                }
            }
        }

        ctx.popCurrentElement();

        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);
        for (const b of builders) ret.add(b.getStructure());

        return ret.getSelection();
    };
}

export interface IncludeConnectedParams {
    query: StructureQuery,
    bondTest?: QueryFn<boolean>,
    layerCount: number,
    wholeResidues: boolean,
    fixedPoint: boolean
}

export function includeConnected({ query, layerCount, wholeResidues, bondTest, fixedPoint }: IncludeConnectedParams): StructureQuery {
    const lc = Math.max(layerCount, 0);
    return function query_includeConnected(ctx) {
        const builder = StructureSelection.UniqueBuilder(ctx.inputStructure);
        const src = query(ctx);
        ctx.pushCurrentBond();
        ctx.atomicBond.setTestFn(bondTest);
        StructureSelection.forEach(src, (s, sI) => {
            let incl = s;

            if (fixedPoint) {
                while (true) {
                    const prevCount = incl.elementCount;
                    incl = includeConnectedStep(ctx, wholeResidues, incl);
                    if (incl.elementCount === prevCount) break;
                }
            } else {
                for (let i = 0; i < lc; i++) {
                    incl = includeConnectedStep(ctx, wholeResidues, incl);
                }
            }
            builder.add(incl);
            if (sI % 10 === 0) ctx.throwIfTimedOut();
        });
        ctx.popCurrentBond();

        return builder.getSelection();
    };
}

function includeConnectedStep(ctx: QueryContext, wholeResidues: boolean, structure: Structure) {
    const expanded = expandConnected(ctx, structure);
    if (wholeResidues) return getWholeResidues(ctx, ctx.inputStructure, expanded);
    return expanded;
}

function expandConnected(ctx: QueryContext, structure: Structure) {
    const inputStructure = ctx.inputStructure;
    const interBonds = inputStructure.interUnitBonds;
    const builder = new StructureUniqueSubsetBuilder(inputStructure);

    const atomicBond = ctx.atomicBond;

    // Process intra unit bonds
    for (const unit of structure.units) {
        if (unit.kind !== Unit.Kind.Atomic) {
            // add the whole unit
            builder.beginUnit(unit.id);
            for (let i = 0, _i = unit.elements.length; i < _i; i++) {
                builder.addElement(unit.elements[i]);
            }
            builder.commitUnit();
            continue;
        }

        const inputUnitA = inputStructure.unitMap.get(unit.id) as Unit.Atomic;
        const { offset: intraBondOffset, b: intraBondB, edgeProps: { flags, order, key } } = inputUnitA.bonds;

        atomicBond.setStructure(inputStructure);

        // Process intra unit bonds
        atomicBond.a.unit = inputUnitA;
        atomicBond.b.unit = inputUnitA;
        for (let i = 0, _i = unit.elements.length; i < _i; i++) {
            // add the current element
            builder.addToUnit(unit.id, unit.elements[i]);

            const aIndex = SortedArray.indexOf(inputUnitA.elements, unit.elements[i]) as StructureElement.UnitIndex;

            // check intra unit bonds
            for (let lI = intraBondOffset[aIndex], _lI = intraBondOffset[aIndex + 1]; lI < _lI; lI++) {
                const bIndex = intraBondB[lI] as StructureElement.UnitIndex;
                const bElement = inputUnitA.elements[bIndex];

                // Check if the element is already present:
                if (SortedArray.has(unit.elements, bElement) || builder.has(unit.id, bElement)) continue;

                atomicBond.aIndex = aIndex;
                atomicBond.a.element = unit.elements[i];
                atomicBond.bIndex = bIndex;
                atomicBond.b.element = bElement;
                atomicBond.type = flags[lI];
                atomicBond.order = order[lI];
                atomicBond.key = key[lI];

                if (atomicBond.test(ctx, true)) {
                    builder.addToUnit(unit.id, bElement);
                }
            }
        }

        // Process inter unit bonds
        for (const bondedUnit of interBonds.getConnectedUnits(inputUnitA.id)) {
            const currentUnitB = structure.unitMap.get(bondedUnit.unitB);
            const inputUnitB = inputStructure.unitMap.get(bondedUnit.unitB) as Unit.Atomic;

            for (const aI of bondedUnit.connectedIndices) {
                // check if the element is in the expanded structure
                if (!SortedArray.has(unit.elements, inputUnitA.elements[aI])) continue;

                for (const bond of bondedUnit.getEdges(aI)) {
                    const bElement = inputUnitB.elements[bond.indexB];

                    // Check if the element is already present:
                    if ((currentUnitB && SortedArray.has(currentUnitB.elements, bElement)) || builder.has(bondedUnit.unitB, bElement)) continue;

                    atomicBond.a.unit = inputUnitA;
                    atomicBond.aIndex = aI;
                    atomicBond.a.element = inputUnitA.elements[aI];
                    atomicBond.b.unit = inputUnitB;
                    atomicBond.bIndex = bond.indexB;
                    atomicBond.b.element = bElement;
                    atomicBond.type = bond.props.flag;
                    atomicBond.order = bond.props.order;
                    atomicBond.key = bond.props.key;

                    if (atomicBond.test(ctx, true)) {
                        builder.addToUnit(bondedUnit.unitB, bElement);
                    }
                }
            }
        }
    }

    return builder.getStructure();
}

export interface SurroundingLigandsParams {
    query: StructureQuery,
    radius: number,
    includeWater: boolean
}

/**
 * Includes expanded surrounding ligands based on radius from the source, struct_conn entries & pdbx_molecule entries.
 */
export function surroundingLigands({ query, radius, includeWater }: SurroundingLigandsParams): StructureQuery {
    return function query_surroundingLigands(ctx) {

        const inner = StructureSelection.unionStructure(query(ctx));
        const surroundings = getWholeResidues(ctx, ctx.inputStructure, getIncludeSurroundings(ctx, ctx.inputStructure, inner, { radius }));

        const prd = getPrdAsymIdx(ctx.inputStructure);
        const graph = getStructConnInfo(ctx.inputStructure);

        const l = StructureElement.Location.create(surroundings);

        const includedPrdChains = new Map<string, string[]>();

        const componentResidues = new ResidueSet({ checkOperator: true });

        for (const unit of surroundings.units) {
            if (unit.kind !== Unit.Kind.Atomic) continue;

            l.unit = unit;

            const { elements } = unit;
            const chainsIt = Segmentation.transientSegments(unit.model.atomicHierarchy.chainAtomSegments, elements);
            const residuesIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements);

            while (chainsIt.hasNext) {
                const chainSegment = chainsIt.move();
                l.element = elements[chainSegment.start];

                const asym_id = StructureProperties.chain.label_asym_id(l);
                const op_name = StructureProperties.unit.operator_name(l);

                // check for PRD molecules
                if (prd.has(asym_id)) {
                    if (includedPrdChains.has(asym_id)) {
                        arraySetAdd(includedPrdChains.get(asym_id)!, op_name);
                    } else {
                        includedPrdChains.set(asym_id, [op_name]);
                    }
                    continue;
                }

                const entityType = StructureProperties.entity.type(l);

                // test entity and chain
                if (entityType === 'water' || entityType === 'polymer') continue;

                residuesIt.setSegment(chainSegment);
                while (residuesIt.hasNext) {
                    const residueSegment = residuesIt.move();
                    l.element = elements[residueSegment.start];
                    graph.addComponent(ResidueSet.getEntryFromLocation(l), componentResidues);
                }
            }

            ctx.throwIfTimedOut();
        }

        // assemble the core structure

        const builder = ctx.inputStructure.subsetBuilder(true);
        for (const unit of ctx.inputStructure.units) {
            if (unit.kind !== Unit.Kind.Atomic) continue;

            l.unit = unit;
            const { elements } = unit;
            const chainsIt = Segmentation.transientSegments(unit.model.atomicHierarchy.chainAtomSegments, elements);
            const residuesIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements);

            builder.beginUnit(unit.id);

            while (chainsIt.hasNext) {
                const chainSegment = chainsIt.move();
                l.element = elements[chainSegment.start];

                const asym_id = StructureProperties.chain.label_asym_id(l);
                const op_name = StructureProperties.unit.operator_name(l);

                if (includedPrdChains.has(asym_id) && includedPrdChains.get(asym_id)!.indexOf(op_name) >= 0) {
                    builder.addElementRange(elements, chainSegment.start, chainSegment.end);
                    continue;
                }

                if (!componentResidues.hasLabelAsymId(asym_id)) {
                    continue;
                }

                residuesIt.setSegment(chainSegment);
                while (residuesIt.hasNext) {
                    const residueSegment = residuesIt.move();
                    l.element = elements[residueSegment.start];

                    if (!componentResidues.has(l)) continue;
                    builder.addElementRange(elements, residueSegment.start, residueSegment.end);
                }
            }
            builder.commitUnit();

            ctx.throwIfTimedOut();
        }

        const components = structureUnion(ctx.inputStructure, [builder.getStructure(), inner]);

        // add water
        if (includeWater) {
            const finalBuilder = new StructureUniqueSubsetBuilder(ctx.inputStructure);
            const lookup = ctx.inputStructure.lookup3d;
            for (const unit of components.units) {
                const c = unit.conformation;
                const elements = unit.elements;
                for (let i = 0, _i = elements.length; i < _i; i++) {
                    const e = elements[i];
                    lookup.findIntoBuilderIf(c.x(e), c.y(e), c.z(e), radius, finalBuilder, testIsWater);
                    finalBuilder.addToUnit(unit.id, e);
                }

                ctx.throwIfTimedOut();
            }

            return StructureSelection.Sequence(ctx.inputStructure, [finalBuilder.getStructure()]);
        } else {
            return StructureSelection.Sequence(ctx.inputStructure, [components]);
        }
    };
}

const _entity_type = StructureProperties.entity.type;
function testIsWater(l: StructureElement.Location) {
    return _entity_type(l) === 'water';
}

function getPrdAsymIdx(structure: Structure) {
    const model = structure.models[0];
    const ids = new Set<string>();
    if (!MmcifFormat.is(model.sourceData)) return ids;
    const { _rowCount, asym_id } = model.sourceData.data.db.pdbx_molecule;
    for (let i = 0; i < _rowCount; i++) {
        ids.add(asym_id.value(i));
    }
    return ids;
}

function getStructConnInfo(structure: Structure) {
    const model = structure.models[0];
    const graph = new StructConnGraph();

    if (!MmcifFormat.is(model.sourceData)) return graph;

    const struct_conn = model.sourceData.data.db.struct_conn;
    const { conn_type_id } = struct_conn;
    const { ptnr1_label_asym_id, ptnr1_label_comp_id, ptnr1_label_seq_id, ptnr1_symmetry, pdbx_ptnr1_label_alt_id, pdbx_ptnr1_PDB_ins_code } = struct_conn;
    const { ptnr2_label_asym_id, ptnr2_label_comp_id, ptnr2_label_seq_id, ptnr2_symmetry, pdbx_ptnr2_label_alt_id, pdbx_ptnr2_PDB_ins_code } = struct_conn;

    for (let i = 0; i < struct_conn._rowCount; i++) {
        const bondType = conn_type_id.value(i);
        if (bondType !== 'covale' && bondType !== 'metalc') continue;

        const a: ResidueSetEntry = {
            label_asym_id: ptnr1_label_asym_id.value(i),
            label_comp_id: ptnr1_label_comp_id.value(i),
            label_seq_id: ptnr1_label_seq_id.value(i),
            label_alt_id: pdbx_ptnr1_label_alt_id.value(i),
            ins_code: pdbx_ptnr1_PDB_ins_code.value(i),
            operator_name: ptnr1_symmetry.value(i) ?? '1_555'
        };

        const b: ResidueSetEntry = {
            label_asym_id: ptnr2_label_asym_id.value(i),
            label_comp_id: ptnr2_label_comp_id.value(i),
            label_seq_id: ptnr2_label_seq_id.value(i),
            label_alt_id: pdbx_ptnr2_label_alt_id.value(i),
            ins_code: pdbx_ptnr2_PDB_ins_code.value(i),
            operator_name: ptnr2_symmetry.value(i) ?? '1_555'
        };

        graph.addEdge(a, b);
    }

    return graph;
}

class StructConnGraph {
    vertices = new Map<string, ResidueSetEntry>();
    edges = new Map<string, string[]>();

    private addVertex(e: ResidueSetEntry, label: string) {
        if (this.vertices.has(label)) return;
        this.vertices.set(label, e);
        this.edges.set(label, []);
    }

    addEdge(a: ResidueSetEntry, b: ResidueSetEntry) {
        const al = ResidueSet.getLabel(a);
        const bl = ResidueSet.getLabel(b);
        this.addVertex(a, al);
        this.addVertex(b, bl);
        arraySetAdd(this.edges.get(al)!, bl);
        arraySetAdd(this.edges.get(bl)!, al);
    }

    addComponent(start: ResidueSetEntry, set: ResidueSet) {
        const startLabel = ResidueSet.getLabel(start);

        if (!this.vertices.has(startLabel)) {
            set.add(start);
            return;
        }

        const visited = new Set<string>();
        const added = new Set<string>();
        const stack = [startLabel];

        added.add(startLabel);
        set.add(start);

        while (stack.length > 0) {
            const a = stack.pop()!;
            visited.add(a);

            const u = this.vertices.get(a)!;

            for (const b of this.edges.get(a)!) {
                if (visited.has(b)) continue;
                stack.push(b);

                if (added.has(b)) continue;
                added.add(b);

                const v = this.vertices.get(b)!;
                if (u.operator_name === v.operator_name) {
                    set.add({ ...v, operator_name: start.operator_name });
                } else {
                    set.add(v);
                }

            }
        }
    }
}

// TODO: unionBy (skip this one?), cluster