/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation } from '../../../../mol-data/int';
import { MinimizeRmsd } from '../../../../mol-math/linear-algebra/3d/minimize-rmsd';
import { BestDatabaseSequenceMapping } from '../../../../mol-model-props/sequence/best-database-mapping';
import { ElementIndex } from '../../model/indexing';
import { Structure } from '../structure';
import { Unit } from '../unit';

export interface AlignmentResult {
    transform: MinimizeRmsd.Result,
    pivot: number,
    other: number
}

export function alignAndSuperposeWithBestDatabaseMapping(structures: Structure[]): AlignmentResult[] {
    const indexMap = new Map<string, IndexEntry>();

    for (let i = 0; i < structures.length; i++) {
        buildIndex(structures[i], indexMap, i);
    }

    const index = Array.from(indexMap.values());

    // TODO: support non-first structure pivots
    const pairs = findPairs(structures.length, index);

    const ret: AlignmentResult[] = [];
    for (const p of pairs) {
        const [a, b] = getPositionTables(index, p.i, p.j, p.count);
        const transform = MinimizeRmsd.compute({ a, b });
        ret.push({ transform, pivot: p.i, other: p.j });
    }

    return ret;
}

function getPositionTables(index: IndexEntry[], pivot: number, other: number, N: number) {
    const xs = MinimizeRmsd.Positions.empty(N);
    const ys = MinimizeRmsd.Positions.empty(N);

    let o = 0;
    for (const { pivots } of index) {
        const a = pivots[pivot];
        const b = pivots[other];
        if (!a || !b) continue;

        const l = Math.min(a[2] - a[1], b[2] - b[1]);

        for (let i = 0; i < l; i++) {
            let eI = (a[1] + i) as ElementIndex;
            xs.x[o] = a[0].conformation.x(eI);
            xs.y[o] = a[0].conformation.y(eI);
            xs.z[o] = a[0].conformation.z(eI);

            eI = (b[1] + i) as ElementIndex;
            ys.x[o] = b[0].conformation.x(eI);
            ys.y[o] = b[0].conformation.y(eI);
            ys.z[o] = b[0].conformation.z(eI);
            o++;
        }
    }

    return [xs, ys];
}

function findPairs(N: number, index: IndexEntry[]) {
    const pairwiseCounts: number[][] = [];
    for (let i = 0; i < N; i++) {
        pairwiseCounts[i] = [];
        for (let j = 0; j < N; j++) pairwiseCounts[i][j] = 0;
    }

    for (const { pivots } of index) {
        for (let i = 0; i < N; i++) {
            if (!pivots[i]) continue;

            const lI = pivots[i]![2] - pivots[i]![1];

            for (let j = i + 1; j < N; j++) {
                if (!pivots[j]) continue;

                const lJ = pivots[j]![2] - pivots[j]![1];
                pairwiseCounts[i][j] = pairwiseCounts[i][j] + Math.min(lI, lJ);
            }
        }
    }

    const ret: { i: number, j: number, count: number }[] = [];

    for (let j = 1; j < N; j++) {
        ret[j - 1] = { i: 0, j, count: pairwiseCounts[0][j] };
    }

    // TODO: support non-first structure pivots
    // for (let i = 0; i < N - 1; i++) {
    //     let max = 0, maxJ = i;
    //     for (let j = i + 1; j < N; j++) {
    //         if (pairwiseCounts[i][j] > max) {
    //             maxJ = j;
    //             max = pairwiseCounts[i][j];
    //         }
    //     }

    //     ret[i] = { i, j: maxJ, count: max };
    // }

    return ret;
}

interface IndexEntry {
    key: string,
    pivots: { [i: number]: [unit: Unit.Atomic, start: ElementIndex, end: ElementIndex] | undefined }
}

function buildIndex(structure: Structure, index: Map<string, IndexEntry>, sI: number) {
    for (const unit of structure.units) {
        if (unit.kind !== Unit.Kind.Atomic) continue;

        const { elements, model } = unit;
        const { offsets: residueOffset } = model.atomicHierarchy.residueAtomSegments;

        const map = BestDatabaseSequenceMapping.Provider.get(model).value;
        if (!map) return;

        const { dbName, accession, num } = map;

        const chainsIt = Segmentation.transientSegments(unit.model.atomicHierarchy.chainAtomSegments, elements);
        const residuesIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements);

        while (chainsIt.hasNext) {
            const chainSegment = chainsIt.move();
            residuesIt.setSegment(chainSegment);
            while (residuesIt.hasNext) {
                const residueSegment = residuesIt.move();
                const eI = elements[residueSegment.start];
                const rI = residueOffset[eI];

                if (!dbName[rI]) continue;

                const key = `${dbName[rI]}-${accession[rI]}-${num[rI]}`;

                if (!index.has(key)) {
                    index.set(key, { key, pivots: { [sI]: [unit, eI, elements[residueSegment.end]] } });
                } else {
                    const entry = index.get(key)!;

                    if (!entry.pivots[sI]) {
                        entry.pivots[sI] = [unit, eI, elements[residueSegment.end]];
                    }
                }
            }
        }
    }
}