/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation, SortedArray } from '../../../../../mol-data/int';
import { IntAdjacencyGraph } from '../../../../../mol-math/graph';
import { BondType } from '../../../model/types';
import { StructureElement } from '../../../structure';
import Unit from '../../unit';
import { IntraUnitBonds } from '../bonds/data';
import { sortArray } from '../../../../../mol-data/util';
import { Column } from '../../../../../mol-data/db';
import { arraySetAdd, arraySetRemove } from '../../../../../mol-util/array';

export function computeRings(unit: Unit.Atomic) {
    const size = largestResidue(unit);
    const state = State(unit, size);

    const residuesIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
    while (residuesIt.hasNext) {
        const seg = residuesIt.move();
        processResidue(state, seg.start, seg.end);
    }

    return state.rings;
}

const enum Constants {
    MaxDepth = 4
}

interface State {
    startVertex: number,
    endVertex: number,
    count: number,
    visited: Int32Array,
    queue: Int32Array,
    color: Int32Array,
    pred: Int32Array,

    left: Int32Array,
    right: Int32Array,

    currentColor: number,
    currentAltLoc: string,
    hasAltLoc: boolean,

    rings: SortedArray<StructureElement.UnitIndex>[],
    currentRings: SortedArray<StructureElement.UnitIndex>[],
    bonds: IntraUnitBonds,
    unit: Unit.Atomic,
    altLoc: Column<string>
}

function State(unit: Unit.Atomic, capacity: number): State {
    return {
        startVertex: 0,
        endVertex: 0,
        count: 0,
        visited: new Int32Array(capacity),
        queue: new Int32Array(capacity),
        pred: new Int32Array(capacity),
        left: new Int32Array(Constants.MaxDepth),
        right: new Int32Array(Constants.MaxDepth),
        color: new Int32Array(capacity),
        currentColor: 0,
        currentAltLoc: '',
        hasAltLoc: false,
        rings: [],
        currentRings: [],
        unit,
        bonds: unit.bonds,
        altLoc: unit.model.atomicHierarchy.atoms.label_alt_id
    };
}

function resetState(state: State) {
    state.count = state.endVertex - state.startVertex;
    const { visited, pred, color } = state;
    for (let i = 0; i < state.count; i++) {
        visited[i] = -1;
        pred[i] = -1;
        color[i] = 0;
    }
    state.currentColor = 0;
    state.currentAltLoc = '';
    state.hasAltLoc = false;
}

function largestResidue(unit: Unit.Atomic) {
    const residuesIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
    let size = 0;
    while (residuesIt.hasNext) {
        const seg = residuesIt.move();
        size = Math.max(size, seg.end - seg.start);
    }
    return size;
}

function processResidue(state: State, start: number, end: number) {
    const { visited } = state;
    state.startVertex = start;
    state.endVertex = end;

    // no two atom rings
    if (state.endVertex - state.startVertex < 3) return;

    state.currentRings = [];

    const { elements } = state.unit;
    const altLocs: string[] = [];
    for (let i = state.startVertex; i < state.endVertex; i++) {
        const altLoc = state.altLoc.value(elements[i]);
        arraySetAdd(altLocs, altLoc);
    }
    arraySetRemove(altLocs, '');

    if (altLocs.length === 0) {
        resetState(state);
        for (let i = 0; i < state.count; i++) {
            if (visited[i] >= 0) continue;
            findRings(state, i);
        }
    } else {
        for (let aI = 0; aI < altLocs.length; aI++) {
            resetState(state);
            state.hasAltLoc = true;
            state.currentAltLoc = altLocs[aI];
            for (let i = 0; i < state.count; i++) {
                if (visited[i] >= 0) continue;
                const altLoc = state.altLoc.value(elements[state.startVertex + i]);
                if (altLoc && altLoc !== state.currentAltLoc) {
                    continue;
                }
                findRings(state, i);
            }
        }
    }

    for (let i = 0, _i = state.currentRings.length; i < _i; i++) {
        state.rings.push(state.currentRings[i]);
    }
}

function addRing(state: State, a: number, b: number) {
    // only "monotonous" rings
    if (b < a) {
        return;
    }

    const { pred, color, left, right } = state;
    const nc = ++state.currentColor;

    let current = a;

    for (let t = 0; t < Constants.MaxDepth; t++) {
        color[current] = nc;
        current = pred[current];
        if (current < 0) break;
    }

    let leftOffset = 0, rightOffset = 0;

    let found = false, target = 0;
    current = b;
    for (let t = 0; t < Constants.MaxDepth; t++) {
        if (color[current] === nc) {
            target = current;
            found = true;
            break;
        }
        right[rightOffset++] = current;
        current = pred[current];
        if (current < 0) break;
    }
    if (!found) {
        return;
    }

    current = a;
    for (let t = 0; t < Constants.MaxDepth; t++) {
        left[leftOffset++] = current;
        if (target === current) break;
        current = pred[current];
        if (current < 0) break;
    }

    const len = leftOffset + rightOffset;
    // rings must have at least three elements
    if (len < 3) {
        return;
    }

    const ring = new Int32Array(len);
    let ringOffset = 0;
    for (let t = 0; t < leftOffset; t++) ring[ringOffset++] = state.startVertex + left[t];
    for (let t = rightOffset - 1; t >= 0; t--) ring[ringOffset++] = state.startVertex + right[t];

    sortArray(ring);

    if (state.hasAltLoc) {
        // we need to check if the ring was already added because alt locs are present.

        for (let rI = 0, _rI = state.currentRings.length; rI < _rI; rI++) {
            const r = state.currentRings[rI];
            if (ring[0] !== r[0]) continue;
            if (ring.length !== r.length) continue;

            let areSame = true;
            for (let aI = 0, _aI = ring.length; aI < _aI; aI++) {
                if (ring[aI] !== r[aI]) {
                    areSame = false;
                    break;
                }
            }
            if (areSame) {
                return;
            }
        }
    }

    state.currentRings.push(SortedArray.ofSortedArray(ring));
}

function findRings(state: State, from: number) {
    const { bonds, startVertex, endVertex, visited, queue, pred } = state;
    const { elements } = state.unit;
    const { b: neighbor, edgeProps: { flags: bondFlags }, offset } = bonds;
    visited[from] = 1;
    queue[0] = from;
    let head = 0, size = 1;

    while (head < size) {
        const top = queue[head++];
        const a = startVertex + top;
        const start = offset[a], end = offset[a + 1];

        for (let i = start; i < end; i++) {
            const b = neighbor[i];
            if (b < startVertex || b >= endVertex || !BondType.isCovalent(bondFlags[i])) continue;

            if (state.hasAltLoc) {
                const altLoc = state.altLoc.value(elements[b]);
                if (altLoc && state.currentAltLoc !== altLoc) {
                    continue;
                }
            }

            const other = b - startVertex;

            if (visited[other] > 0) {
                if (pred[other] !== top && pred[top] !== other) {
                    addRing(state, top, other);
                }
                continue;
            }

            visited[other] = 1;
            queue[size++] = other;
            pred[other] = top;
        }
    }
}

export function getFingerprint(elements: string[]) {
    const len = elements.length;
    const reversed: string[] = new Array(len);

    for (let i = 0; i < len; i++) reversed[i] = elements[len - i - 1];

    const rotNormal = getMinimalRotation(elements);
    const rotReversed = getMinimalRotation(reversed);

    let isNormalSmaller = false;

    for (let i = 0; i < len; i++) {
        const u = elements[(i + rotNormal) % len], v = reversed[(i + rotReversed) % len];
        if (u !== v) {
            isNormalSmaller = u < v;
            break;
        }
    }

    if (isNormalSmaller) return buildFinderprint(elements, rotNormal);
    return buildFinderprint(reversed, rotReversed);
}

function getMinimalRotation(elements: string[]) {
    // adapted from http://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation

    const len = elements.length;
    const f = new Int32Array(len * 2);
    for (let i = 0; i < f.length; i++) f[i] = -1;

    let u = '', v = '', k = 0;

    for (let j = 1; j < f.length; j++) {
        let i = f[j - k - 1];
        while (i !== -1) {
            u = elements[j % len]; v = elements[(k + i + 1) % len];
            if (u === v) break;
            if (u < v) k = j - i - 1;
            i = f[i];
        }

        if (i === -1) {
            u = elements[j % len]; v = elements[(k + i + 1) % len];
            if (u !== v) {
                if (u < v) k = j;
                f[j - k] = -1;
            } else f[j - k] = i + 1;
        } else f[j - k] = i + 1;
    }

    return k;
}

function buildFinderprint(elements: string[], offset: number) {
    const len = elements.length;
    const ret: string[] = [];
    let i;
    for (i = 0; i < len - 1; i++) {
        ret.push(elements[(i + offset) % len]);
        ret.push('-');
    }
    ret.push(elements[(i + offset) % len]);
    return ret.join('');
}

type RingIndex = import('../rings').UnitRings.Index
type RingComponentIndex = import('../rings').UnitRings.ComponentIndex

export function createIndex(rings: ArrayLike<SortedArray<StructureElement.UnitIndex>>, aromaticRings: ReadonlyArray<RingIndex>) {
    const elementRingIndices: Map<StructureElement.UnitIndex, RingIndex[]> = new Map();
    const elementAromaticRingIndices: Map<StructureElement.UnitIndex, RingIndex[]> = new Map();

    // for each ring atom, assign all rings that it is present in
    for (let rI = 0 as RingIndex, _rI = rings.length; rI < _rI; rI++) {
        const r = rings[rI];
        for (let i = 0, _i = r.length; i < _i; i++) {
            const e = r[i];
            if (elementRingIndices.has(e)) elementRingIndices.get(e)!.push(rI);
            else elementRingIndices.set(e, [rI]);
        }
    }

    // for each ring atom, assign all aromatic rings that it is present in
    for (let aI = 0, _aI = aromaticRings.length; aI < _aI; aI++) {
        const rI = aromaticRings[aI];
        const r = rings[rI];
        for (let i = 0, _i = r.length; i < _i; i++) {
            const e = r[i];
            if (elementAromaticRingIndices.has(e)) elementAromaticRingIndices.get(e)!.push(rI);
            else elementAromaticRingIndices.set(e, [rI]);
        }
    }

    // create a graph where vertices are rings, edge if two rings share at least one atom
    const graph = new IntAdjacencyGraph.UniqueEdgeBuilder(rings.length);
    for (let rI = 0 as RingIndex, _rI = rings.length; rI < _rI; rI++) {
        const r = rings[rI];

        for (let i = 0, _i = r.length; i < _i; i++) {
            const e = r[i];

            const containedRings = elementRingIndices.get(e)!;

            if (containedRings.length === 1) continue;

            for (let j = 0, _j = containedRings.length; j < _j; j++) {
                const rJ = containedRings[j];
                if (rI >= rJ) continue;
                graph.addEdge(rI, rJ);
            }
        }
    }

    const components = IntAdjacencyGraph.connectedComponents(graph.getGraph());

    const ringComponentIndex = components.componentIndex as any as RingComponentIndex[];
    const ringComponents: RingIndex[][] = [];
    for (let i = 0; i < components.componentCount; i++) ringComponents[i] = [];

    for (let rI = 0 as RingIndex, _rI = rings.length; rI < _rI; rI++) {
        ringComponents[ringComponentIndex[rI]].push(rI);
    }

    return { elementRingIndices, elementAromaticRingIndices, ringComponentIndex, ringComponents };
}