/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, StructureElement, Unit } from '../../mol-model/structure';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { PluginContext } from '../../mol-plugin/context';
import { StateBuilder, StateObjectCell, StateSelection, StateTransform } from '../../mol-state';
import { StructureComponentRef } from '../manager/structure/hierarchy-state';
import { EmptyLoci, isEmptyLoci, Loci } from '../../mol-model/loci';
import { Wiggle } from '../../mol-theme/wiggle';
import { OrderedSet } from '../../mol-data/int';

type WiggleEachReprCallback = (update: StateBuilder.Root, repr: StateObjectCell<PluginStateObject.Molecule.Structure.Representation3D, StateTransform<typeof StateTransforms.Representation.StructureRepresentation3D>>, wiggle?: StateObjectCell<any, StateTransform<typeof StateTransforms.Representation.WiggleStructureRepresentation3DFromBundle>>) => Promise<void>
const WiggleManagerTag = 'wiggle-controls';

export async function setStructureWiggle(plugin: PluginContext, components: StructureComponentRef[], value: number, lociGetter: (structure: Structure) => Promise<StructureElement.Loci | EmptyLoci>, types?: string[]) {
    await eachRepr(plugin, components, async (update, repr, wiggleCell) => {
        if (types && types.length > 0 && !types.includes(repr.params!.values.type.name)) return;

        const structure = repr.obj!.data.sourceData;
        // always use the root structure to get the loci so the wiggle
        // stays applicable as long as the root structure does not change
        const loci = await lociGetter(structure.root);
        if (Loci.isEmpty(loci) || isEmptyLoci(loci)) return;

        const layer = {
            bundle: StructureElement.Bundle.fromLoci(loci),
            value,
        };

        if (wiggleCell) {
            const bundleLayers = [...wiggleCell.params!.values.layers, layer];
            const filtered = getFilteredBundle(bundleLayers, structure);
            update.to(wiggleCell).update(Wiggle.toBundle(filtered));
        } else {
            const filtered = getFilteredBundle([layer], structure);
            update.to(repr.transform.ref)
                .apply(StateTransforms.Representation.WiggleStructureRepresentation3DFromBundle, Wiggle.toBundle(filtered), { tags: WiggleManagerTag });
        }
    });
}

export async function clearStructureWiggle(plugin: PluginContext, components: StructureComponentRef[], types?: string[]) {
    await eachRepr(plugin, components, async (update, repr, wiggleCell) => {
        if (types && types.length > 0 && !types.includes(repr.params!.values.type.name)) return;
        if (wiggleCell) {
            update.delete(wiggleCell.transform.ref);
        }
    });
}

async function eachRepr(plugin: PluginContext, components: StructureComponentRef[], callback: WiggleEachReprCallback) {
    const state = plugin.state.data;
    const update = state.build();
    for (const c of components) {
        for (const r of c.representations) {
            const wiggle = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.WiggleStructureRepresentation3DFromBundle, r.cell.transform.ref).withTag(WiggleManagerTag));
            await callback(update, r.cell, wiggle[0]);
        }
    }

    return update.commit({ doNotUpdateCurrent: true });
}

/** filter wiggle layers for given structure */
function getFilteredBundle(layers: Wiggle.BundleLayer[], structure: Structure) {
    const wiggle = Wiggle.ofBundle(layers, structure.root);
    const merged = Wiggle.merge(wiggle);
    return Wiggle.filter(merged, structure) as Wiggle<StructureElement.Loci>;
}

function getUncertaintyValue(unit: Unit, element: number): number {
    if (Unit.isAtomic(unit)) {
        return unit.model.atomicConformation.B_iso_or_equiv.value(element);
    } else if (Unit.isSpheres(unit)) {
        return unit.model.coarseConformation.spheres.rmsf[element];
    }
    return 0;
}

/** Compute min/max uncertainty (B-factor or RMSF) across all units in a structure */
function getUncertaintyRange(structure: Structure): { min: number, max: number } {
    let min = Infinity;
    let max = -Infinity;
    for (const unit of structure.units) {
        const elements = unit.elements;
        for (let j = 0, jl = elements.length; j < jl; j++) {
            const v = getUncertaintyValue(unit, elements[j]);
            if (v < min) min = v;
            if (v > max) max = v;
        }
    }
    if (!isFinite(min)) min = 0;
    if (!isFinite(max)) max = 0;
    return { min, max };
}

/**
 * Set per-group wiggle based on B-factor/RMSF uncertainty data.
 * Values are normalized to [0, 1] within each structure's min-max range.
 * @param scale - maximum wiggle value (default 1.0, corresponds to Angstroms when combined with wiggleAmplitude)
 */
export async function setStructureWiggleFromUncertainty(plugin: PluginContext, components: StructureComponentRef[], scale: number = 1, types?: string[]) {
    await eachRepr(plugin, components, async (update, repr, wiggleCell) => {
        if (types && types.length > 0 && !types.includes(repr.params!.values.type.name)) return;

        const structure = repr.obj!.data.sourceData;
        const root = structure.root;
        const { min, max } = getUncertaintyRange(root);
        const range = max - min;
        if (range <= 0) return;

        // Group elements by discretized uncertainty bucket (256 levels for Uint8 texture)
        const buckets = new Map<number, { unit: Unit, indices: number[] }[]>();

        for (const unit of root.units) {
            const elements = unit.elements;
            const unitBuckets = new Map<number, number[]>();

            for (let j = 0, jl = elements.length; j < jl; j++) {
                const v = getUncertaintyValue(unit, elements[j]);
                const normalized = (v - min) / range;
                const bucket = Math.min(255, Math.round(normalized * 255));
                if (!unitBuckets.has(bucket)) unitBuckets.set(bucket, []);
                unitBuckets.get(bucket)!.push(j);
            }

            for (const [bucket, indices] of unitBuckets) {
                if (!buckets.has(bucket)) buckets.set(bucket, []);
                buckets.get(bucket)!.push({ unit, indices });
            }
        }

        // Create one layer per bucket
        const bundleLayers: Wiggle.BundleLayer[] = [];
        for (const [bucket, unitIndices] of buckets) {
            const value = (bucket / 255) * scale;
            const elements: StructureElement.Loci['elements'][0][] = [];
            for (const { unit, indices } of unitIndices) {
                elements.push({
                    unit,
                    indices: OrderedSet.ofSortedArray(new Int32Array(indices) as any as StructureElement.UnitIndex[]),
                });
            }
            const loci = StructureElement.Loci(root, elements);
            if (!StructureElement.Loci.isEmpty(loci)) {
                bundleLayers.push({
                    bundle: StructureElement.Bundle.fromLoci(loci),
                    value,
                });
            }
        }

        if (bundleLayers.length === 0) return;

        const filtered = getFilteredBundle(bundleLayers, structure);
        if (wiggleCell) {
            update.to(wiggleCell).update(Wiggle.toBundle(filtered));
        } else {
            update.to(repr.transform.ref)
                .apply(StateTransforms.Representation.WiggleStructureRepresentation3DFromBundle, Wiggle.toBundle(filtered), { tags: WiggleManagerTag });
        }
    });
}
