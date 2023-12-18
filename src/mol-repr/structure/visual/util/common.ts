/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit, Structure, ElementIndex, StructureElement, ResidueIndex } from '../../../../mol-model/structure';
import { Mat4, Vec3 } from '../../../../mol-math/linear-algebra';
import { TransformData, createTransform } from '../../../../mol-geo/geometry/transform-data';
import { OrderedSet, SortedArray } from '../../../../mol-data/int';
import { EmptyLoci, Loci } from '../../../../mol-model/loci';
import { AtomicNumbers } from '../../../../mol-model/structure/model/properties/atomic';
import { fillSerial } from '../../../../mol-util/array';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { AssignableArrayLike } from '../../../../mol-util/type-helpers';
import { getBoundary } from '../../../../mol-math/geometry/boundary';
import { Box3D, Sphere3D } from '../../../../mol-math/geometry';
import { SizeTheme } from '../../../../mol-theme/size';
import { hasPolarNeighbour } from '../../../../mol-model-props/computed/chemistry/functional-group';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const m4toArray = Mat4.toArray;

/** Return a Loci for the elements of a whole residue the elementIndex belongs to. */
export function getResidueLoci(structure: Structure, unit: Unit.Atomic, elementIndex: ElementIndex): Loci {
    const { elements, model } = unit;
    if (OrderedSet.indexOf(elements, elementIndex) !== -1) {
        const { index, offsets } = model.atomicHierarchy.residueAtomSegments;
        const rI = index[elementIndex];
        const _indices: number[] = [];
        for (let i = offsets[rI], il = offsets[rI + 1]; i < il; ++i) {
            const unitIndex = OrderedSet.indexOf(elements, i);
            if (unitIndex !== -1) _indices.push(unitIndex);
        }
        const indices = OrderedSet.ofSortedArray<StructureElement.UnitIndex>(SortedArray.ofSortedArray(_indices));
        return StructureElement.Loci(structure, [{ unit, indices }]);
    }
    return EmptyLoci;
}

/**
 * Return a Loci for the elements of a whole residue the elementIndex belongs to but
 * restrict to elements that have the same label_alt_id or none
 */
export function getAltResidueLoci(structure: Structure, unit: Unit.Atomic, elementIndex: ElementIndex) {
    const { elements, model } = unit;
    const { label_alt_id } = model.atomicHierarchy.atoms;
    const elementAltId = label_alt_id.value(elementIndex);
    if (OrderedSet.indexOf(elements, elementIndex) !== -1) {
        const { index } = model.atomicHierarchy.residueAtomSegments;
        const rI = index[elementIndex];
        return getAltResidueLociFromId(structure, unit, rI, elementAltId);
    }
    return StructureElement.Loci(structure, []);
}

export function getAltResidueLociFromId(structure: Structure, unit: Unit.Atomic, residueIndex: ResidueIndex, elementAltId: string) {
    const { elements, model } = unit;
    const { label_alt_id } = model.atomicHierarchy.atoms;
    const { offsets } = model.atomicHierarchy.residueAtomSegments;

    const _indices: number[] = [];
    for (let i = offsets[residueIndex], il = offsets[residueIndex + 1]; i < il; ++i) {
        const unitIndex = OrderedSet.indexOf(elements, i);
        if (unitIndex !== -1) {
            const altId = label_alt_id.value(i);
            if (elementAltId === altId || altId === '') {
                _indices.push(unitIndex);
            }
        }
    }
    const indices = OrderedSet.ofSortedArray<StructureElement.UnitIndex>(SortedArray.ofSortedArray(_indices));
    return StructureElement.Loci(structure, [{ unit, indices }]);
}

//

export type StructureGroup = { structure: Structure, group: Unit.SymmetryGroup }

export function createUnitsTransform(structureGroup: StructureGroup, includeParent: boolean, invariantBoundingSphere: Sphere3D, cellSize: number, batchSize: number, transformData?: TransformData) {
    const { child } = structureGroup.structure;
    const units: ReadonlyArray<Unit> = includeParent && child
        ? structureGroup.group.units.filter(u => child.unitMap.has(u.id))
        : structureGroup.group.units;
    const unitCount = units.length;
    const n = unitCount * 16;
    const array = transformData && transformData.aTransform.ref.value.length >= n ? transformData.aTransform.ref.value : new Float32Array(n);
    for (let i = 0; i < unitCount; i++) {
        m4toArray(units[i].conformation.operator.matrix, array, i * 16);
    }
    return createTransform(array, unitCount, invariantBoundingSphere, cellSize, batchSize, transformData);
}

export const UnitKindInfo = {
    'atomic': {},
    'spheres': {},
    'gaussians': {},
};
export type UnitKind = keyof typeof UnitKindInfo
export const UnitKindOptions = PD.objectToOptions(UnitKindInfo);

export function includesUnitKind(unitKinds: UnitKind[], unit: Unit) {
    for (let i = 0, il = unitKinds.length; i < il; ++i) {
        if (Unit.isAtomic(unit) && unitKinds[i] === 'atomic') return true;
        if (Unit.isSpheres(unit) && unitKinds[i] === 'spheres') return true;
        if (Unit.isGaussians(unit) && unitKinds[i] === 'gaussians') return true;
    }
    return false;
}

//

const DefaultMaxCells = 500_000_000;

export function getVolumeSliceInfo(box: Box3D, resolution: number, maxCells = DefaultMaxCells) {
    const size = Box3D.size(Vec3(), box);
    Vec3.ceil(size, size);
    size.sort((a, b) => b - a); // descending
    const maxAreaCells = Math.floor(Math.cbrt(maxCells) * Math.cbrt(maxCells));
    const area = size[0] * size[1];
    const areaCells = Math.ceil(area / (resolution * resolution));
    return { area, areaCells, maxAreaCells };
}

/**
 * Guard against overly high resolution for the given box size.
 * Internally it uses the largest 2d slice of the box to determine the
 * maximum resolution to account for the 2d texture layout on the GPU.
 */
export function ensureReasonableResolution<T>(box: Box3D, props: { resolution: number } & T, maxCells = DefaultMaxCells) {
    const { area, areaCells, maxAreaCells } = getVolumeSliceInfo(box, props.resolution, maxCells);
    const resolution = areaCells > maxAreaCells ? Math.sqrt(area / maxAreaCells) : props.resolution;
    return { ...props, resolution };
}

export function getConformation(unit: Unit) {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return unit.model.atomicConformation;
        case Unit.Kind.Spheres: return unit.model.coarseConformation.spheres;
        case Unit.Kind.Gaussians: return unit.model.coarseConformation.gaussians;
    }
}

export const CommonSurfaceParams = {
    ignoreHydrogens: PD.Boolean(false, { description: 'Whether or not to include hydrogen atoms in the surface calculation.' }),
    ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
    traceOnly: PD.Boolean(false, { description: 'Whether or not to only use trace atoms in the surface calculation.' }),
    includeParent: PD.Boolean(false, { description: 'Include elements of the parent structure in surface calculation to get a surface patch of the current structure.' }),
};
export const DefaultCommonSurfaceProps = PD.getDefaultValues(CommonSurfaceParams);
export type CommonSurfaceProps = typeof DefaultCommonSurfaceProps

const v = Vec3();
function squaredDistance(x: number, y: number, z: number, center: Vec3) {
    return Vec3.squaredDistance(Vec3.set(v, x, y, z), center);
}

/** marks `indices` for filtering/ignoring in `id` when not in `elements` */
function filterUnitId(id: AssignableArrayLike<number>, elements: SortedArray, indices: SortedArray) {
    let start = 0;
    const end = elements.length;
    for (let i = 0, il = indices.length; i < il; ++i) {
        const idx = SortedArray.indexOfInRange(elements, indices[i], start, end);
        if (idx === -1) {
            id[i] = -2;
        } else {
            id[i] = idx;
            start = idx;
        }
    }
}

export function getUnitConformationAndRadius(structure: Structure, unit: Unit, sizeTheme: SizeTheme<any>, props: CommonSurfaceProps) {
    const { ignoreHydrogens, ignoreHydrogensVariant, traceOnly, includeParent } = props;
    const rootUnit = includeParent ? structure.root.unitMap.get(unit.id) : unit;
    const differentRoot = includeParent && rootUnit !== unit;

    const { x, y, z } = getConformation(rootUnit);
    const { elements } = rootUnit;
    const { center, radius: sphereRadius } = unit.boundary.sphere;
    const extraRadius = (4 + 1.5) * 2; // TODO should be twice (the max vdW/sphere radius plus the probe radius)
    const radiusSq = (sphereRadius + extraRadius) * (sphereRadius + extraRadius);

    let indices: SortedArray<ElementIndex>;
    let id: AssignableArrayLike<number>;

    if (ignoreHydrogens || traceOnly || differentRoot) {
        const _indices: number[] = [];
        const _id: number[] = [];
        for (let i = 0, il = elements.length; i < il; ++i) {
            const eI = elements[i];
            if (ignoreHydrogens && isHydrogen(structure, rootUnit, eI, ignoreHydrogensVariant)) continue;
            if (traceOnly && !isTrace(rootUnit, eI)) continue;
            if (differentRoot && squaredDistance(x[eI], y[eI], z[eI], center) > radiusSq) continue;

            _indices.push(eI);
            _id.push(i);
        }
        indices = SortedArray.ofSortedArray(_indices);
        id = _id;
    } else {
        indices = elements;
        id = fillSerial(new Int32Array(indices.length));
    }

    if (includeParent && rootUnit !== unit) {
        filterUnitId(id, unit.elements, indices);
    }

    const position = { indices, x, y, z, id };
    const boundary = differentRoot ? getBoundary(position) : unit.boundary;

    const l = StructureElement.Location.create(structure, rootUnit);
    const radius = (index: number) => {
        l.element = index as ElementIndex;
        return sizeTheme.size(l);
    };

    return { position, boundary, radius };
}

export function getStructureConformationAndRadius(structure: Structure, sizeTheme: SizeTheme<any>, props: CommonSurfaceProps) {
    const { ignoreHydrogens, ignoreHydrogensVariant, traceOnly, includeParent } = props;
    const differentRoot = includeParent && !!structure.parent;
    const l = StructureElement.Location.create(structure.root);

    const { center, radius: sphereRadius } = structure.boundary.sphere;
    const extraRadius = (4 + 1.5) * 2; // TODO should be twice (the max vdW/sphere radius plus the probe radius)
    const radiusSq = (sphereRadius + extraRadius) * (sphereRadius + extraRadius);

    let xs: ArrayLike<number>;
    let ys: ArrayLike<number>;
    let zs: ArrayLike<number>;
    let rs: ArrayLike<number>;
    let id: AssignableArrayLike<number>;
    let indices: OrderedSet<number>;

    if (ignoreHydrogens || traceOnly || differentRoot) {
        const { getSerialIndex } = structure.serialMapping;
        const units = differentRoot ? structure.root.units : structure.units;

        const _xs: number[] = [];
        const _ys: number[] = [];
        const _zs: number[] = [];
        const _rs: number[] = [];
        const _id: number[] = [];
        for (let i = 0, il = units.length; i < il; ++i) {
            const unit = units[i];
            const { elements, conformation: c } = unit;
            const childUnit = structure.unitMap.get(unit.id);

            l.unit = unit;
            for (let j = 0, jl = elements.length; j < jl; ++j) {
                const eI = elements[j];
                if (ignoreHydrogens && isHydrogen(structure, unit, eI, ignoreHydrogensVariant)) continue;
                if (traceOnly && !isTrace(unit, eI)) continue;

                const _x = c.x(eI), _y = c.y(eI), _z = c.z(eI);
                if (differentRoot && squaredDistance(_x, _y, _z, center) > radiusSq) continue;

                _xs.push(_x);
                _ys.push(_y);
                _zs.push(_z);
                l.element = eI;
                _rs.push(sizeTheme.size(l));

                if (differentRoot) {
                    const idx = childUnit ? SortedArray.indexOf(childUnit.elements, eI) : -1;
                    if (idx === -1) {
                        _id.push(-2); // mark for filtering/ignoring when not in `elements`
                    } else {
                        _id.push(getSerialIndex(childUnit, eI));
                    }
                } else {
                    _id.push(getSerialIndex(unit, eI));
                }
            }
        }
        xs = _xs, ys = _ys, zs = _zs, rs = _rs;
        id = _id;
        indices = OrderedSet.ofRange(0, id.length);
    } else {
        const { elementCount } = structure;
        const _xs = new Float32Array(elementCount);
        const _ys = new Float32Array(elementCount);
        const _zs = new Float32Array(elementCount);
        const _rs = new Float32Array(elementCount);
        for (let i = 0, m = 0, il = structure.units.length; i < il; ++i) {
            const unit = structure.units[i];
            const { elements, conformation: c } = unit;
            l.unit = unit;
            for (let j = 0, jl = elements.length; j < jl; ++j) {
                const eI = elements[j];

                const mj = m + j;
                _xs[mj] = c.x(eI);
                _ys[mj] = c.y(eI);
                _zs[mj] = c.z(eI);
                l.element = eI;
                _rs[mj] = sizeTheme.size(l);
            }
            m += elements.length;
        }
        xs = _xs, ys = _ys, zs = _zs, rs = _rs;
        id = fillSerial(new Uint32Array(elementCount));
        indices = OrderedSet.ofRange(0, id.length);
    }

    const position = { indices, x: xs, y: ys, z: zs, id };
    const boundary = differentRoot ? getBoundary(position) : structure.boundary;
    const radius = (index: number) => rs[index];

    return { position, boundary, radius };
}

const _H = AtomicNumbers['H'];
export function isHydrogen(structure: Structure, unit: Unit, element: ElementIndex, variant: 'all' | 'non-polar' | 'polar') {
    if (Unit.isCoarse(unit)) return false;
    if (unit.model.atomicHierarchy.derived.atom.atomicNumber[element] !== _H) return false;
    if (variant === 'all') return true;
    const polar = hasPolarNeighbour(structure, unit, SortedArray.indexOf(unit.elements, element) as StructureElement.UnitIndex);
    if (polar && variant === 'polar') return true;
    if (!polar && variant === 'non-polar') return true;
    return false;
}
export function isH(atomicNumber: ArrayLike<number>, element: ElementIndex) {
    return atomicNumber[element] === _H;
}

export function isTrace(unit: Unit, element: ElementIndex) {
    if (Unit.isCoarse(unit)) return true;
    const atomId = unit.model.atomicHierarchy.atoms.label_atom_id.value(element);
    if (atomId === 'CA' || atomId === 'BB' || atomId === 'P') return true;
    return false;
}
