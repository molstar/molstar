/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SphericalBasisOrder } from '../../extensions/alpha-orbitals/spherical-functions';
import { CreateOrbitalRepresentation3D, CreateOrbitalVolume, StaticBasisAndOrbitals } from '../../extensions/alpha-orbitals/transforms';
import { createPluginAsync, DefaultPluginSpec } from '../../mol-plugin';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectSelector, StateTransformer } from '../../mol-state';
import { Color } from '../../mol-util/color';
import { ColorNames } from '../../mol-util/color/names';
import { ParamDefinition } from '../../mol-util/param-definition';
import { mountControls } from './controls';
import { DemoMoleculeSDF, DemoOrbitals } from './example-data';
import { BehaviorSubject } from 'rxjs';
import { debounceTime, skip } from 'rxjs/operators';
import './index.html';
import { Basis, AlphaOrbital } from '../../extensions/alpha-orbitals/data-model';
import { canComputeAlphaOrbitalsOnGPU } from '../../extensions/alpha-orbitals/gpu/compute';
import { PluginCommands } from '../../mol-plugin/commands';
require('mol-plugin-ui/skin/light.scss');

interface DemoInput {
    moleculeSdf: string,
    basis: Basis,
    order: SphericalBasisOrder,
    orbitals: AlphaOrbital[]
}

interface Params {
    orbitalIndex: number,
    isoValue: number,
    gpuSurface: boolean
}

export class AlphaOrbitalsExample {
    plugin: PluginContext;

    async init(target: string | HTMLElement) {
        this.plugin = await createPluginAsync(typeof target === 'string' ? document.getElementById(target)! : target, {
            ...DefaultPluginSpec,
            layout: {
                initial: {
                    isExpanded: false,
                    showControls: false
                },
                controls: { left: 'none', right: 'none', top: 'none', bottom: 'none' },
            },
            config: [
                [PluginConfig.Viewport.ShowExpand, false],
                [PluginConfig.Viewport.ShowControls, false],
                [PluginConfig.Viewport.ShowSelectionMode, false],
                [PluginConfig.Viewport.ShowAnimation, false],
            ]
        });

        this.plugin.managers.interactivity.setProps({ granularity: 'element' });

        if (!canComputeAlphaOrbitalsOnGPU(this.plugin.canvas3d?.webgl)) {
            PluginCommands.Toast.Show(this.plugin, {
                title: 'Error',
                message: `Browser/device does not support required WebGL extension (OES_texture_float).`
            });
            return;
        }

        this.load({
            moleculeSdf: DemoMoleculeSDF,
            ...DemoOrbitals
        });

        mountControls(this, document.getElementById('controls')!);
    }

    readonly params = new BehaviorSubject<ParamDefinition.For<Params>>({} as any);
    readonly state = new BehaviorSubject<Params>({ orbitalIndex: 32, isoValue: 1, gpuSurface: false });

    private volume?: StateObjectSelector<PluginStateObject.Volume.Data, typeof CreateOrbitalVolume>;
    private positive?: StateObjectSelector<PluginStateObject.Volume.Representation3D, typeof CreateOrbitalRepresentation3D>;
    private negative?: StateObjectSelector<PluginStateObject.Volume.Representation3D, typeof CreateOrbitalRepresentation3D>;
    private currentParams: Params = { ...this.state.value };

    private async setIndex() {
        if (!this.volume?.isOk) return;
        const state = this.state.value;
        await this.plugin.build().to(this.volume).update(CreateOrbitalVolume, () => ({ index: state.orbitalIndex })).commit();
        if (this.currentParams.gpuSurface !== this.state.value.gpuSurface) {
            await this.setIsovalue();
        }
        this.currentParams = this.state.value;
    }

    private setIsovalue() {
        this.currentParams = this.state.value;
        const update = this.plugin.build();
        update.to(this.positive!).update(this.volumeParams('positive', ColorNames.blue));
        update.to(this.negative!).update(this.volumeParams('negative', ColorNames.red));
        return update.commit();
    }

    private volumeParams(kind: 'positive' | 'negative', color: Color): StateTransformer.Params<typeof CreateOrbitalRepresentation3D> {
        return {
            alpha: 0.85,
            color,
            directVolume: this.state.value.gpuSurface,
            kind,
            relativeIsovalue: this.state.value.isoValue,
            pickable: false,
            xrayShaded: true
        };
    }

    async load(input: DemoInput) {
        await this.plugin.clear();

        const data = await this.plugin.builders.data.rawData({ data: input.moleculeSdf }, { state: { isGhost: true } });
        const trajectory = await this.plugin.builders.structure.parseTrajectory(data, 'mol');
        const model = await this.plugin.builders.structure.createModel(trajectory);
        const structure = await this.plugin.builders.structure.createStructure(model);

        const all = await this.plugin.builders.structure.tryCreateComponentStatic(structure, 'all');
        if (all) await this.plugin.builders.structure.representation.addRepresentation(all, { type: 'ball-and-stick', color: 'element-symbol', colorParams: { carbonColor: { name: 'element-symbol', params: {} } } });

        const state = this.state.value;

        this.volume = await this.plugin.build().toRoot()
            .apply(StaticBasisAndOrbitals, { basis: input.basis, order: input.order, orbitals: input.orbitals })
            .apply(CreateOrbitalVolume, { index: state.orbitalIndex })
            .commit();

        if (!this.volume.isOk) {
            this.volume = void 0;
            return;
        }

        const repr = this.plugin.build().to(this.volume);

        this.positive = repr.apply(CreateOrbitalRepresentation3D, this.volumeParams('positive', ColorNames.blue)).selector;
        this.negative = repr.apply(CreateOrbitalRepresentation3D, this.volumeParams('negative', ColorNames.red)).selector;

        await repr.commit();

        this.params.next({
            orbitalIndex: ParamDefinition.Numeric(this.currentParams.orbitalIndex, { min: 0, max: input.orbitals.length - 1 }, { immediateUpdate: true, isEssential: true }),
            isoValue: ParamDefinition.Numeric(this.currentParams.isoValue, { min: 0.5, max: 3, step: 0.1 }, { immediateUpdate: true, isEssential: false }),
            gpuSurface: ParamDefinition.Boolean(this.currentParams.gpuSurface)
        });

        this.state.pipe(skip(1), debounceTime(1000 / 24)).subscribe(async params => {
            if (params.orbitalIndex !== this.currentParams.orbitalIndex) {
                this.setIndex();
            } else if (params.isoValue !== this.currentParams.isoValue || params.gpuSurface !== this.currentParams.gpuSurface) {
                this.setIsovalue();
            }
        });
    }
}

(window as any).AlphaOrbitalsExample = new AlphaOrbitalsExample();