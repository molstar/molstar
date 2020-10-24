/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Basis, computeIsocontourValues } from '../../extensions/alpha-orbitals/cubes';
import { SphericalBasisOrder } from '../../extensions/alpha-orbitals/orbitals';
import { createPluginAsync, DefaultPluginSpec } from '../../mol-plugin';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { PluginContext } from '../../mol-plugin/context';
import { ColorNames } from '../../mol-util/color/names';
import { DemoMoleculeSDF, DemoOrbitals } from './example-data';
import './index.html';
import { CreateOrbitalVolume, StaticBasisAndOrbitals } from './transforms';
require('mol-plugin-ui/skin/light.scss');

interface DemoInput {
    moleculeSdf: string,
    basis: Basis,
    order: SphericalBasisOrder,
    orbitals: {
        energy: number,
        alpha: number[]
    }[]
}

class AlphaOrbitalsExample {
    plugin: PluginContext;

    async init(target: string | HTMLElement) {
        this.plugin = await createPluginAsync(typeof target === 'string' ? document.getElementById(target)! : target, {
            ...DefaultPluginSpec,
            layout: {
                initial: {
                    isExpanded: false,
                    showControls: false
                },
                controls: { left: 'none', right: 'none', top: 'none', bottom: 'none' }
            }
        });

        this.plugin.managers.interactivity.setProps({ granularity: 'element' });

        this.load({
            moleculeSdf: DemoMoleculeSDF,
            ...DemoOrbitals
        });
    }

    async load(input: DemoInput) {
        await this.plugin.clear();

        const data = await this.plugin.builders.data.rawData({ data: input.moleculeSdf }, { state: { isGhost: true } });
        const trajectory = await this.plugin.builders.structure.parseTrajectory(data, 'mol');
        const model = await this.plugin.builders.structure.createModel(trajectory);
        const structure = await this.plugin.builders.structure.createStructure(model);

        const all = await this.plugin.builders.structure.tryCreateComponentStatic(structure, 'all');
        if (all) await this.plugin.builders.structure.representation.addRepresentation(all, { type: 'ball-and-stick', color: 'element-symbol', colorParams: { carbonColor: { name: 'element-symbol', params: { } } } });

        const volumeRef = await this.plugin.build().toRoot()
            .apply(StaticBasisAndOrbitals, { basis: input.basis, order: input.order, orbitals: input.orbitals })
            .apply(CreateOrbitalVolume, { index: 44 })
            .commit();

        if (!volumeRef.isOk) return;

        // TODO: the isovalues are being computed twice. Need to add more flexible support to Volume object
        //       for controlling them
        const { negative, positive } = computeIsocontourValues(volumeRef.data!.grid.cells.data as any, 0.85);

        const repr = this.plugin.build().to(volumeRef);


        if (positive !== void 0) {
            repr.apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(this.plugin, volumeRef.data!, {
                type: 'isosurface',
                typeParams: { isoValue: { kind: 'absolute', absoluteValue: positive } },
                color: 'uniform',
                colorParams: { value: ColorNames.blue }
            }));
        }

        if (negative !== void 0) {
            repr.apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(this.plugin, volumeRef.data!, {
                type: 'isosurface',
                typeParams: { isoValue: { kind: 'absolute', absoluteValue: negative } },
                color: 'uniform',
                colorParams: { value: ColorNames.red }
            }));
        }

        await repr.commit();
    }
}

(window as any).AlphaOrbitalsExample = new AlphaOrbitalsExample();