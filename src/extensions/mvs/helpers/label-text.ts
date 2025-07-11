/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Sphere3D } from '../../../mol-math/geometry';
import { BoundaryHelper } from '../../../mol-math/geometry/boundary-helper';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { ElementIndex, Model, Structure, StructureElement, StructureProperties, Unit } from '../../../mol-model/structure';
import { arrayExtend } from '../../../mol-util/array';
import { AtomRanges } from './atom-ranges';
import { IndicesAndSortings } from './indexing';
import { MVSAnnotationRow } from './schemas';
import { getAtomRangesForRows } from './selections';
import { isDefined } from './utils';


/** Properties describing position, size, etc. of a text in 3D */
export interface TextProps {
    /** Anchor point for the text (i.e. the center of the text will appear in front of `center`) */
    center: Vec3,
    /** Depth of the text wrt anchor point (i.e. the text will appear in distance `radius` in front of the anchor point) */
    depth: number,
    /** Relative text size */
    scale: number,
    /** Index of the first atom within structure, to which this text is bound (for coloring and similar purposes) */
    group: number,
}

const tmpVec = Vec3();
const tmpArray: number[] = [];
const boundaryHelper = new BoundaryHelper('98');
const outAtoms: ElementIndex[] = [];
const outFirstAtomIndex: { value?: number } = {};

/** Helper for caching atom ranges qualifying to a group of annotation rows, per `Unit`. */
class AtomRangesCache {
    private readonly cache: { [key: string]: AtomRanges } = {};
    private readonly hasOperators: boolean;

    constructor(private readonly rows: MVSAnnotationRow[]) {
        this.hasOperators = rows.some(row => isDefined(row.instance_id));
    }

    get(unit: Unit): AtomRanges {
        const instanceId = unit.conformation.operator.instanceId;
        const key = this.hasOperators ? `${unit.model.id}:${instanceId}` : unit.model.id;
        return this.cache[key] ??= getAtomRangesForRows(this.rows, unit.model, instanceId, IndicesAndSortings.get(unit.model));
    }
}

/** Return `TextProps` (position, size, etc.) for a text that is to be bound to a substructure of `structure` defined by union of `rows`.
 * Derives `center` and `depth` from the boundary sphere of the substructure, `scale` from the number of heavy atoms in the substructure. */
export function textPropsForSelection(structure: Structure, sizeFunction: (location: StructureElement.Location) => number, rows: MVSAnnotationRow[], onlyInModel?: Model): TextProps | undefined {
    const loc = StructureElement.Location.create(structure);
    const { units } = structure;
    const { type_symbol } = StructureProperties.atom;
    tmpArray.length = 0;
    let includedAtoms = 0;
    let includedHeavyAtoms = 0;
    let group: number | undefined = undefined;
    let atomSize: number | undefined = undefined;
    const atomRangesCache = new AtomRangesCache(rows);
    for (let iUnit = 0, nUnits = units.length; iUnit < nUnits; iUnit++) {
        const unit = units[iUnit];
        if (onlyInModel && unit.model.id !== onlyInModel.id) continue;
        const ranges = atomRangesCache.get(unit);
        loc.unit = unit;
        AtomRanges.selectAtomsInRanges(unit.elements, ranges, outAtoms, outFirstAtomIndex);
        for (const atom of outAtoms) {
            loc.element = atom;
            unit.conformation.position(atom, tmpVec);
            arrayExtend(tmpArray, tmpVec);
            group ??= structure.serialMapping.cumulativeUnitElementCount[iUnit] + outFirstAtomIndex.value!;
            atomSize ??= sizeFunction(loc);
            includedAtoms++;
            if (type_symbol(loc) !== 'H') includedHeavyAtoms++;
        }
    }
    if (includedAtoms > 0) {
        const { center, radius } = (includedAtoms > 1) ? boundarySphere(tmpArray) : { center: Vec3.fromArray(Vec3(), tmpArray, 0), radius: 1.1 * atomSize! };
        const scale = (includedHeavyAtoms || includedAtoms) ** (1 / 3);
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
