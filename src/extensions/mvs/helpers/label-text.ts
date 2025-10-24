/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Sphere3D } from '../../../mol-math/geometry';
import { BoundaryHelper } from '../../../mol-math/geometry/boundary-helper';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { ElementIndex, Model, Structure, StructureElement, StructureProperties, Unit } from '../../../mol-model/structure';
import { getPhysicalRadius } from '../../../mol-theme/size/physical';
import { arrayExtend } from '../../../mol-util/array';
import { ElementRanges } from './element-ranges';
import { IndicesAndSortings } from './indexing';
import { MVSAnnotationRow } from './schemas';
import { getAtomRangesForRows, getGaussianRangesForRows, getSphereRangesForRows } from './selections';
import { isDefined } from './utils';


/** Properties describing position, size, etc. of a text in 3D */
export interface TextProps {
    /** Anchor point for the text (i.e. the center of the text will appear in front of `center`) */
    center: Vec3,
    /** Depth of the text wrt anchor point (i.e. the text will appear in distance `radius` in front of the anchor point) */
    depth: number,
    /** Relative text size */
    scale: number,
    /** Index of the first element within structure, to which this text is bound (for coloring and similar purposes) */
    group: number,
}

const tmpVec = Vec3();
const tmpArray: number[] = [];
const boundaryHelper = new BoundaryHelper('98');
const outElements: ElementIndex[] = [];
const outFirstElementIndex: { value?: number } = {};

/** Helper for caching element ranges qualifying to a group of annotation rows, per `Unit`. */
class ElementRangesCache {
    private readonly cache: { [key: string]: ElementRanges } = {};
    private readonly hasOperators: boolean;

    constructor(private readonly rows: MVSAnnotationRow[]) {
        this.hasOperators = rows.some(row => isDefined(row.instance_id));
    }

    get(unit: Unit): ElementRanges {
        const instanceId = unit.conformation.operator.instanceId;
        const key = `${unit.model.id}:${unit.kind}:${this.hasOperators ? instanceId : '*'}`;
        return this.cache[key] ??= this.compute(unit);
    }
    private compute(unit: Unit): ElementRanges {
        const instanceId = unit.conformation.operator.instanceId;
        const indices = IndicesAndSortings.get(unit.model);
        switch (unit.kind) {
            case Unit.Kind.Atomic:
                return getAtomRangesForRows(this.rows, unit.model, instanceId, indices);
            case Unit.Kind.Spheres:
                return getSphereRangesForRows(this.rows, unit.model, instanceId, indices);
            case Unit.Kind.Gaussians:
                return getGaussianRangesForRows(this.rows, unit.model, instanceId, indices);
        }
    }
}

/** Approximate number of heavy atoms per protein residue (I got 7.55 from 2e2n) */
const AVG_ATOMS_PER_RESIDUE = 8;

/** Return `TextProps` (position, size, etc.) for a text that is to be bound to a substructure of `structure` defined by union of `rows`.
 * Derives `center` and `depth` from the boundary sphere of the substructure, `scale` from the number of heavy atoms in the substructure. */
export function textPropsForSelection(structure: Structure, rows: MVSAnnotationRow[], onlyInModel?: Model): TextProps | undefined {
    const loc = StructureElement.Location.create(structure);
    const { units } = structure;
    const { type_symbol } = StructureProperties.atom;
    tmpArray.length = 0;
    let includedElements = 0;
    let includedHeavyAtoms = 0;
    let group: number | undefined = undefined;
    /** Used for `depth` in case the selection has only 1 element (hence bounding sphere radius is 0) */
    let singularRadius: number | undefined = undefined;
    const elementRangesCache = new ElementRangesCache(rows);
    for (let iUnit = 0, nUnits = units.length; iUnit < nUnits; iUnit++) {
        const unit = units[iUnit];
        if (onlyInModel && unit.model.id !== onlyInModel.id) continue;
        const coarseElements = unit.kind === Unit.Kind.Spheres ? unit.model.coarseHierarchy.spheres : unit.kind === Unit.Kind.Gaussians ? unit.model.coarseHierarchy.gaussians : undefined;
        const ranges = elementRangesCache.get(unit);
        ElementRanges.selectElementsInRanges(unit.elements, ranges, outElements, outFirstElementIndex);
        loc.unit = unit;
        for (const iElem of outElements) {
            loc.element = iElem;
            arrayExtend(tmpArray, unit.conformation.position(iElem, tmpVec));
            group ??= structure.serialMapping.cumulativeUnitElementCount[iUnit] + outFirstElementIndex.value!;
            singularRadius ??= getPhysicalRadius(unit, iElem) * 1.2;
            if (coarseElements) {
                // coarse
                const nResidues = coarseElements.seq_id_end.value(iElem) - coarseElements.seq_id_begin.value(iElem) + 1;
                includedHeavyAtoms += nResidues * AVG_ATOMS_PER_RESIDUE;
            } else {
                // atomic
                if (type_symbol(loc) !== 'H') includedHeavyAtoms++;
            }
            includedElements++;
        }
    }
    if (includedElements > 0) {
        const { center, radius } = (includedElements > 1) ? boundarySphere(tmpArray) : { center: Vec3.fromArray(Vec3(), tmpArray, 0), radius: singularRadius! };
        const scale = (includedHeavyAtoms || includedElements) ** (1 / 3);
        return { center, depth: radius, scale, group: group! };
    }
}

/** Calculate the boundary sphere for a set of points given by their flattened coordinates (`flatCoords.slice(0,3)` is the first point etc.) */
function boundarySphere(flatCoords: readonly number[]): Sphere3D {
    const length = flatCoords.length;
    boundaryHelper.reset();
    for (let offset = 0; offset < length; offset += 3) {
        Vec3.fromArray(tmpVec, flatCoords, offset);
        boundaryHelper.includePosition(tmpVec);
    }
    boundaryHelper.finishedIncludeStep();
    for (let offset = 0; offset < length; offset += 3) {
        Vec3.fromArray(tmpVec, flatCoords, offset);
        boundaryHelper.radiusPosition(tmpVec);
    }
    return boundaryHelper.getSphere();
}
