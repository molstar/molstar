/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, StructureElement } from '../structure';
import { StructureProperties } from '../structure/structure/properties';

/**
 * Build a per-target sub-structure map from a parent structure and the
 * `targetMapping` found on a `ParticleList`.
 *
 * For each target ID → chain IDs entry in `targetMapping`, the corresponding
 * units in `parentStructure` whose `label_asym_id` is listed are collected into
 * a new `Structure`. Targets with no matching units are skipped.
 */
export function buildTargetStructuresFromMapping(
    parentStructure: Structure,
    targetMapping: ReadonlyMap<number, ReadonlyArray<string>>
): Map<number, Structure> {
    const result = new Map<number, Structure>();

    for (const [targetId, asymIds] of targetMapping) {
        const asymSet = new Set(asymIds);
        const builder = Structure.Builder({ label: parentStructure.label });

        for (const unit of parentStructure.units) {
            const loc = StructureElement.Location.create(parentStructure, unit, unit.elements[0]);
            const asym = StructureProperties.chain.label_asym_id(loc);
            if (asymSet.has(asym)) {
                builder.addUnit(unit.kind, unit.model, unit.conformation.operator, unit.elements, unit.traits);
            }
        }

        const sub = builder.getStructure();
        if (sub.elementCount > 0) {
            result.set(targetId, sub);
        }
    }

    return result;
}
