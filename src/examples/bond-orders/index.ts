/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { BondOrders, BondOrdersTrajectoryPreset } from '../../extensions/bond-orders';
import { Structure, Unit } from '../../mol-model/structure';
import { hasIntraBondOrderFromTable } from '../../mol-model/structure/model/properties/atomic/bonds';
import { WaterNames } from '../../mol-model/structure/model/types';
import { BuiltInTrajectoryFormat } from '../../mol-plugin-state/formats/trajectory';
import { createPluginUI } from '../../mol-plugin-ui';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import { PluginSpec } from '../../mol-plugin/spec';
import { Asset } from '../../mol-util/assets';
import '../../mol-plugin-ui/skin/light.scss';
import './index.html';

type LoadParams = { url: string, format?: BuiltInTrajectoryFormat, isBinary?: boolean }

function logPerceivedDoubles(structure: Structure) {
    const { type_symbol, label_comp_id } = structure.models[0].atomicHierarchy.atoms;
    let doubles = 0;
    for (const unit of structure.units) {
        if (!Unit.isAtomic(unit)) continue;
        const { offset, b, edgeProps: { order } } = unit.bonds;
        const { elements } = unit;
        for (let u = 0, ul = elements.length; u < ul; u++) {
            const compId = label_comp_id.value(elements[u]);
            if (WaterNames.has(compId) || hasIntraBondOrderFromTable(compId)) continue;
            for (let i = offset[u], il = offset[u + 1]; i < il; i++) {
                if (u >= b[i]) continue;
                if (order[i] === 2 && type_symbol.value(elements[u]) === 'C' && type_symbol.value(elements[b[i]]) === 'C') {
                    doubles++;
                }
            }
        }
    }
    console.log(`Bond-orders example: ${doubles} perceived intra-residue C=C (placeholder)`);
}

class BondOrdersExample {
    async init(target: string | HTMLElement) {
        const spec = DefaultPluginUISpec();
        const plugin = await createPluginUI({
            target: typeof target === 'string' ? document.getElementById(target)! : target,
            render: renderReact18,
            spec: {
                ...spec,
                layout: {
                    initial: {
                        isExpanded: false,
                        showControls: true
                    }
                },
                components: {
                    remoteState: 'none'
                },
                behaviors: [
                    ...spec.behaviors,
                    PluginSpec.Behavior(BondOrders)
                ]
            }
        });

        const load = async ({ url, format = 'pdb', isBinary = false }: LoadParams) => {
            await plugin.clear();
            const data = await plugin.builders.data.download({ url: Asset.Url(url), isBinary }, { state: { isGhost: true } });
            const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
            const hierarchy = await plugin.builders.structure.hierarchy.applyPreset(trajectory, BondOrdersTrajectoryPreset);
            const structure = hierarchy?.structure.obj?.data;
            if (structure) logPerceivedDoubles(structure);
        };

        await load({ url: 'https://files.rcsb.org/download/1CBS.pdb', format: 'pdb' });

        (window as any).loadBondOrdersExample = load;
        return plugin;
    }
}

(window as any).BondOrdersExample = new BondOrdersExample();
