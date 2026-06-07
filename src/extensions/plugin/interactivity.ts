/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement } from '../../mol-model/structure';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { MolScriptBuilder } from '../../mol-script/language/builder';
import { Expression } from '../../mol-script/language/expression';
import { CameraFocusLociOptions } from '../../mol-plugin-state/manager/camera';
import { PluginContext } from '../../mol-plugin/context';

export interface StructureInteractivityOptions {
    expression?: (queryBuilder: typeof MolScriptBuilder) => Expression,
    elements?: StructureElement.Schema,
    action: 'highlight' | 'select' | 'focus' | ('highlight' | 'select' | 'focus')[],
    applyGranularity?: boolean,
    filterStructure?: (structure: Structure) => boolean,
    focusOptions?: Partial<CameraFocusLociOptions>
}

/**
 * Triggers structure element selection or highlighting based on the provided
 * MolScript expression or StructureElement schema. Focus action will only apply to the
 * first structure that matches the criteria.
 *
 * If neither `expression` nor `elements` are provided, all selections/highlights
 * will be cleared based on the specified `action`.
 */
export function applyStructureInteractivity(plugin: PluginContext, { expression, elements, action: action_, applyGranularity = false, filterStructure, focusOptions }: StructureInteractivityOptions) {
    const actions = Array.isArray(action_) ? action_ : [action_];

    if (!expression && !elements) {
        if (actions.includes('select')) {
            plugin.managers.interactivity.lociSelects.deselectAll();
        }
        if (actions.includes('highlight')) {
            plugin.managers.interactivity.lociHighlights.clearHighlights();
        }
        return;
    }

    if (actions.includes('select')) {
        plugin.managers.interactivity.lociSelects.deselectAll();
    }

    const structures = plugin.state.data.selectQ(Q => Q.rootsOfType(PluginStateObject.Molecule.Structure));
    let focused = false;
    for (const s of structures) {
        if (!s.obj?.data) continue;

        if (filterStructure && !filterStructure(s.obj.data)) continue;

        const loci = expression
            ? StructureElement.Loci.fromExpression(s.obj.data, expression)
            : StructureElement.Loci.fromSchema(s.obj.data, elements!);

        for (const action of actions) {
            if (action === 'select') {
                plugin.managers.interactivity.lociSelects.select({ loci }, applyGranularity);
            } else if (action === 'highlight') {
                plugin.managers.interactivity.lociHighlights.highlight({ loci }, applyGranularity);
            } else if (action === 'focus' && !StructureElement.Loci.isEmpty(loci) && !focused) {
                plugin.managers.camera.focusLoci(loci, focusOptions);
                focused = true;
                if (actions.length === 1) return; // if only focusing, focus the first matching structure and return immediately
            }
        }
    }
}