/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PresetProvider } from '../preset-provider';
import { PluginStateObject } from '../../objects';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StateObjectRef, StateTransformer } from '../../../mol-state';
import { StateTransforms } from '../../transforms';
import { RootStructureDefinition } from '../../helpers/root-structure';
import { PresetStructureRepresentations, StructureRepresentationPresetProvider } from './representation-preset';
import { PluginContext } from '../../../mol-plugin/context';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { Model, Structure } from '../../../mol-model/structure';
import { getStructureQuality } from '../../../mol-repr/util';
import { OperatorNameColorThemeProvider } from '../../../mol-theme/color/operator-name';
import { PluginConfig } from '../../../mol-plugin/config';
import { CCDFormat } from '../../../mol-model-formats/structure/mmcif';
import { MinimizeRmsd } from '../../../mol-math/linear-algebra/3d/minimize-rmsd';
import { SetUtils } from '../../../mol-util/set';

export interface TrajectoryHierarchyPresetProvider<P = any, S = {}> extends PresetProvider<PluginStateObject.Molecule.Trajectory, P, S> { }
export function TrajectoryHierarchyPresetProvider<P, S>(preset: TrajectoryHierarchyPresetProvider<P, S>) { return preset; }
export namespace TrajectoryHierarchyPresetProvider {
    export type Params<P extends TrajectoryHierarchyPresetProvider> = P extends TrajectoryHierarchyPresetProvider<infer T> ? T : never;
    export type State<P extends TrajectoryHierarchyPresetProvider> = P extends TrajectoryHierarchyPresetProvider<infer _, infer S> ? S : never;

    export const CommonParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
        modelProperties: PD.Optional(PD.Group(StateTransformer.getParamDefinition(StateTransforms.Model.CustomModelProperties, void 0, plugin))),
        structureProperties: PD.Optional(PD.Group(StateTransformer.getParamDefinition(StateTransforms.Model.CustomStructureProperties, void 0, plugin))),
        representationPreset: PD.Optional(PD.Text<keyof PresetStructureRepresentations>('auto' as const))
    });
}

const CommonParams = TrajectoryHierarchyPresetProvider.CommonParams;

const DefaultParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
    model: PD.Optional(PD.Group(StateTransformer.getParamDefinition(StateTransforms.Model.ModelFromTrajectory, a, plugin))),
    showUnitcell: PD.Optional(PD.Boolean(false)),
    structure: PD.Optional(RootStructureDefinition.getParams(void 0, 'assembly').type),
    representationPresetParams: PD.Optional(PD.Group(StructureRepresentationPresetProvider.CommonParams)),
    ...CommonParams(a, plugin)
});

const defaultPreset = TrajectoryHierarchyPresetProvider({
    id: 'preset-trajectory-default',
    display: {
        name: 'Default (Assembly)', group: 'Preset',
        description: 'Shows the first assembly or, if that is unavailable, the first model.'
    },
    isApplicable: o => {
        return true;
    },
    params: DefaultParams,
    async apply(trajectory, params, plugin) {
        const builder = plugin.builders.structure;

        const model = await builder.createModel(trajectory, params.model);
        const modelProperties = await builder.insertModelProperties(model, params.modelProperties);

        const structure = await builder.createStructure(modelProperties || model, params.structure);
        const structureProperties = await builder.insertStructureProperties(structure, params.structureProperties);

        const unitcell = params.showUnitcell === void 0 || !!params.showUnitcell ? await builder.tryCreateUnitcell(modelProperties, undefined, { isHidden: true }) : void 0;
        const representationPreset = params.representationPreset || plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) || PresetStructureRepresentations.auto.id;
        const representation = await plugin.builders.structure.representation.applyPreset(structureProperties, representationPreset, params.representationPresetParams);

        return {
            model,
            modelProperties,
            unitcell,
            structure,
            structureProperties,
            representation
        };
    }
});

const AllModelsParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
    useDefaultIfSingleModel: PD.Optional(PD.Boolean(false)),
    representationPresetParams: PD.Optional(PD.Group(StructureRepresentationPresetProvider.CommonParams)),
    ...CommonParams(a, plugin)
});

const allModels = TrajectoryHierarchyPresetProvider({
    id: 'preset-trajectory-all-models',
    display: {
        name: 'All Models', group: 'Preset',
        description: 'Shows all models; colored by trajectory-index.'
    },
    isApplicable: o => {
        return o.data.frameCount > 1;
    },
    params: AllModelsParams,
    async apply(trajectory, params, plugin) {
        const tr = StateObjectRef.resolveAndCheck(plugin.state.data, trajectory)?.obj?.data;
        if (!tr) return { };

        if (tr.frameCount === 1 && params.useDefaultIfSingleModel) {
            return defaultPreset.apply(trajectory, params as any, plugin);
        }

        const builder = plugin.builders.structure;

        const models = [], structures = [];

        for (let i = 0; i < tr.frameCount; i++) {
            const model = await builder.createModel(trajectory, { modelIndex: i });
            const modelProperties = await builder.insertModelProperties(model, params.modelProperties, { isCollapsed: true });
            const structure = await builder.createStructure(modelProperties || model, { name: 'model', params: {} });
            const structureProperties = await builder.insertStructureProperties(structure, params.structureProperties);

            models.push(model);
            structures.push(structure);

            const quality = structure.obj ? getStructureQuality(structure.obj.data, { elementCountFactor: tr.frameCount }) : 'medium';
            const representationPreset = params.representationPreset || plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) || PresetStructureRepresentations.auto.id;
            await builder.representation.applyPreset(structureProperties, representationPreset, { theme: { globalName: 'trajectory-index' }, quality });
        }

        return { models, structures };
    }
});

const CCDParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
    representationPresetParams: PD.Optional(PD.Group(StructureRepresentationPresetProvider.CommonParams)),
    showOriginalCoordinates: PD.Optional(PD.Boolean(true, { description: `Show original coordinates for 'model' and 'ideal' structure and do not align them.` })),
    shownCoordinateType: PD.Select('ideal', PD.arrayToOptions(['ideal', 'model', 'both'] as const), { description: `What coordinate sets are visible.` }),
    ...CommonParams(a, plugin)
});

const ccd = TrajectoryHierarchyPresetProvider({
    id: 'preset-trajectory-ccd',
    display: {
        name: 'Chemical Component', group: 'Preset',
        description: 'Shows molecules from the Chemical Component Dictionary.'
    },
    isApplicable: o => {
        return CCDFormat.is(o.data.representative.sourceData);
    },
    params: CCDParams,
    async apply(trajectory, params, plugin) {
        const tr = StateObjectRef.resolveAndCheck(plugin.state.data, trajectory)?.obj?.data;
        if (!tr) return {};

        const builder = plugin.builders.structure;

        const idealModel = await builder.createModel(trajectory, { modelIndex: 0 });
        const idealModelProperties = await builder.insertModelProperties(idealModel, params.modelProperties, { isCollapsed: true });

        const idealStructure = await builder.createStructure(idealModelProperties || idealModel, { name: 'model', params: {} });
        const idealStructureProperties = await builder.insertStructureProperties(idealStructure, params.structureProperties);

        const representationPreset = params.representationPreset || PresetStructureRepresentations['chemical-component'].id;
        const representationPresetParams = params.representationPresetParams || {};
        if (representationPresetParams.ignoreHydrogens === undefined) representationPresetParams.ignoreHydrogens = true;

        // degenerate case where either model or ideal coordinates are missing
        if (tr.frameCount !== 2) {
            // 'ideal' references 1st model but it might actually be 'model'
            const coordinateType = Model.CCDCoordinateType.get(idealModel.obj?.data!).coordinateType;
            await builder.representation.applyPreset(idealStructureProperties, representationPreset, { ...representationPresetParams, coordinateType });

            return { models: [idealModel], structures: [idealStructure] };
        }

        const modelModel = await builder.createModel(trajectory, { modelIndex: 1 });
        const modelModelProperties = await builder.insertModelProperties(modelModel, params.modelProperties, { isCollapsed: true });

        const modelStructure = await builder.createStructure(modelModelProperties || modelModel, { name: 'model', params: {} });
        const modelStructureProperties = await builder.insertStructureProperties(modelStructure, params.structureProperties);

        // align ideal and model coordinates
        if (!params.showOriginalCoordinates) {
            const [a, b] = getPositionTables(idealStructure.obj!.data, modelStructure.obj!.data);
            if (!a) {
                plugin.log.warn(`Cannot align ligands whose atom sets are disjoint.`);
            } else {
                const { bTransform, rmsd } = MinimizeRmsd.compute({ a, b });
                await transform(plugin, modelStructure.cell!, bTransform);
                plugin.log.info(`Superposed [model] and [ideal] with RMSD ${rmsd.toFixed(2)}.`);
            }
        }

        await builder.representation.applyPreset(idealStructureProperties, representationPreset, { ...representationPresetParams, coordinateType: CCDFormat.CoordinateType.Ideal, isHidden: params.shownCoordinateType === 'model' });
        await builder.representation.applyPreset(modelStructureProperties, representationPreset, { ...representationPresetParams, coordinateType: CCDFormat.CoordinateType.Model, isHidden: params.shownCoordinateType === 'ideal' });

        return { models: [idealModel, modelModel], structures: [idealStructure, modelStructure] };
    }
});

/** tailored to CCD structures and not generally applicable */
function getPositionTables(s1: Structure, s2: Structure) {
    const m1 = getAtomIdSerialMap(s1);
    const m2 = getAtomIdSerialMap(s2);
    const intersecting = SetUtils.intersection(new Set(m1.keys()), new Set(m2.keys()));

    const ret = [
        MinimizeRmsd.Positions.empty(intersecting.size),
        MinimizeRmsd.Positions.empty(intersecting.size)
    ];
    let o = 0;
    intersecting.forEach(k => {
        ret[0].x[o] = s1.model.atomicConformation.x[m1.get(k)!];
        ret[0].y[o] = s1.model.atomicConformation.y[m1.get(k)!];
        ret[0].z[o] = s1.model.atomicConformation.z[m1.get(k)!];
        ret[1].x[o] = s2.model.atomicConformation.x[m2.get(k)!];
        ret[1].y[o] = s2.model.atomicConformation.y[m2.get(k)!];
        ret[1].z[o] = s2.model.atomicConformation.z[m2.get(k)!];
        o++;
    });

    return ret;
}

/** tailored to CCD structures and not generally applicable */
function getAtomIdSerialMap(structure: Structure) {
    const map = new Map<string, number>();
    const { label_atom_id } = structure.model.atomicHierarchy.atoms;
    for (let i = 0, il = label_atom_id.rowCount; i < il; ++i) {
        const id = label_atom_id.value(i);
        if (!map.has(id)) map.set(id, map.size);
    }
    return map;
}

function transform(plugin: PluginContext, s: StateObjectRef<PluginStateObject.Molecule.Structure>, matrix: Mat4) {
    const b = plugin.state.data.build().to(s)
        .insert(StateTransforms.Model.TransformStructureConformation, { transform: { name: 'matrix', params: { data: matrix, transpose: false } } });
    return plugin.runTask(plugin.state.data.updateTree(b));
}

const CrystalSymmetryParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
    model: PD.Optional(PD.Group(StateTransformer.getParamDefinition(StateTransforms.Model.ModelFromTrajectory, a, plugin))),
    ...CommonParams(a, plugin)
});

async function applyCrystalSymmetry(props: { ijkMin: Vec3, ijkMax: Vec3, theme?: string }, trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory>, params: PD.ValuesFor<ReturnType<typeof CrystalSymmetryParams>>, plugin: PluginContext) {
    const builder = plugin.builders.structure;

    const model = await builder.createModel(trajectory, params.model);
    const modelProperties = await builder.insertModelProperties(model, params.modelProperties);

    const structure = await builder.createStructure(modelProperties || model, {
        name: 'symmetry',
        params: props
    });
    const structureProperties = await builder.insertStructureProperties(structure, params.structureProperties);

    const unitcell = await builder.tryCreateUnitcell(modelProperties, undefined, { isHidden: false });
    const representationPreset = params.representationPreset || plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) || PresetStructureRepresentations.auto.id;
    const representation = await plugin.builders.structure.representation.applyPreset(structureProperties, representationPreset, { theme: { globalName: props.theme } });

    return {
        model,
        modelProperties,
        unitcell,
        structure,
        structureProperties,
        representation
    };
}

const unitcell = TrajectoryHierarchyPresetProvider({
    id: 'preset-trajectory-unitcell',
    display: {
        name: 'Unit Cell', group: 'Preset',
        description: 'Shows the fully populated unit cell.'
    },
    isApplicable: o => {
        return Model.hasCrystalSymmetry(o.data.representative);
    },
    params: CrystalSymmetryParams,
    async apply(trajectory, params, plugin) {
        return await applyCrystalSymmetry({ ijkMin: Vec3.create(0, 0, 0), ijkMax: Vec3.create(0, 0, 0) }, trajectory, params, plugin);
    }
});

const supercell = TrajectoryHierarchyPresetProvider({
    id: 'preset-trajectory-supercell',
    display: {
        name: 'Super Cell', group: 'Preset',
        description: 'Shows the super cell, i.e. the central unit cell and all adjacent unit cells.'
    },
    isApplicable: o => {
        return Model.hasCrystalSymmetry(o.data.representative);
    },
    params: CrystalSymmetryParams,
    async apply(trajectory, params, plugin) {
        return await applyCrystalSymmetry({ ijkMin: Vec3.create(-1, -1, -1), ijkMax: Vec3.create(1, 1, 1), theme: 'operator-hkl' }, trajectory, params, plugin);
    }
});

const CrystalContactsParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
    model: PD.Optional(PD.Group(StateTransformer.getParamDefinition(StateTransforms.Model.ModelFromTrajectory, a, plugin))),
    ...CommonParams(a, plugin)
});

const crystalContacts = TrajectoryHierarchyPresetProvider({
    id: 'preset-trajectory-crystal-contacts',
    display: {
        name: 'Crystal Contacts', group: 'Preset',
        description: 'Showsasymetric unit and chains from neighbours within 5 \u212B, i.e., symmetry mates.'
    },
    isApplicable: o => {
        return Model.hasCrystalSymmetry(o.data.representative);
    },
    params: CrystalContactsParams,
    async apply(trajectory, params, plugin) {
        const builder = plugin.builders.structure;

        const model = await builder.createModel(trajectory, params.model);
        const modelProperties = await builder.insertModelProperties(model, params.modelProperties);

        const structure = await builder.createStructure(modelProperties || model, {
            name: 'symmetry-mates',
            params: { radius: 5 }
        });
        const structureProperties = await builder.insertStructureProperties(structure, params.structureProperties);

        const unitcell = await builder.tryCreateUnitcell(modelProperties, undefined, { isHidden: true });
        const representationPreset = params.representationPreset || plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) || PresetStructureRepresentations.auto.id;
        const representation = await plugin.builders.structure.representation.applyPreset(structureProperties, representationPreset, { theme: { globalName: 'operator-name', carbonColor: 'operator-name', focus: { name: 'element-symbol', params: { carbonColor: { name: 'operator-name', params: OperatorNameColorThemeProvider.defaultValues } } } } });

        return {
            model,
            modelProperties,
            unitcell,
            structure,
            structureProperties,
            representation
        };
    }
});

export const PresetTrajectoryHierarchy = {
    'default': defaultPreset,
    'all-models': allModels,
    ccd,
    unitcell,
    supercell,
    crystalContacts,
};
export type PresetTrajectoryHierarchy = typeof PresetTrajectoryHierarchy;