/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject, PluginStateTransform } from '../../mol-plugin-state/objects';
import { Basis, createSphericalCollocationGrid } from '../../extensions/alpha-orbitals/cubes';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Task } from '../../mol-task';
import { CustomProperties } from '../../mol-model/custom-property';
import { SphericalBasisOrder } from '../../extensions/alpha-orbitals/orbitals';
import { Volume } from '../../mol-model/volume';

export class BasisAndOrbitals extends PluginStateObject.Create<{ basis: Basis, order: SphericalBasisOrder, orbitals: { energy: number, alpha: number[] }[] }>({ name: 'Basis', typeClass: 'Object' }) { }

export const StaticBasisAndOrbitals = PluginStateTransform.BuiltIn({
    name: 'static-basis-and-orbitals',
    display: 'Basis and Orbitals',
    from: PluginStateObject.Root,
    to: BasisAndOrbitals,
    params: {
        label: PD.Text('Orbital Data', { isHidden: true }),
        basis: PD.Value<Basis>(void 0 as any, { isHidden: true }),
        order: PD.Text<SphericalBasisOrder>('gaussian' as SphericalBasisOrder, { isHidden: true }),
        orbitals: PD.Value<{ energy: number, alpha: number[] }[]>([], { isHidden: true })
    },
})({
    apply({ params }) {
        return new BasisAndOrbitals({ basis: params.basis, order: params.order, orbitals: params.orbitals }, { label: params.label });
    }
});

export const CreateOrbitalVolume = PluginStateTransform.BuiltIn({
    name: 'create-orbital-volume',
    display: 'Orbital Volume',
    from: BasisAndOrbitals,
    to: PluginStateObject.Volume.Data,
    params: (a) => {
        if (!a) {
            return { index: PD.Numeric(0) };
        }

        return {
            index: PD.Select(0, a.data.orbitals.map((o, i) => [i, `[${i + 1}] ${o.energy.toFixed(4)}`]))
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Orbital Volume', async ctx => {
            const data = await createSphericalCollocationGrid({
                basis: a.data.basis,
                cutoffThreshold: 0.0075,
                alphaOrbitals: a.data.orbitals[params.index].alpha,
                sphericalOrder: a.data.order,
                boxExpand: 4.5,
                gridSpacing: [
                    [55, 0.6],
                    [40, 0.5],
                    [25, 0.4],
                    [0, 0.33],
                ],
            }).runInContext(ctx);
            const volume: Volume = {
                grid: data.grid,
                sourceData: { name: 'custom grid', kind: 'custom', data },
                customProperties: new CustomProperties(),
                _propertyData: Object.create(null),
            };

            return new PluginStateObject.Volume.Data(volume, { label: 'Orbital Volume' });
        });
    }
});