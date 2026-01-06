/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 */
import { Asset } from '../../mol-util/assets';
import { PluginUIContext } from '../../mol-plugin-ui/context';
import { RibocodeMmcifParseParams, RibocodeMmcifProvider } from './mol-plugin-state/formats/trajectory';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectRef } from '../../mol-state';
import { AlignmentData } from './types';

// Type representing the result of applying a structure preset.
export type PresetResult = Awaited<
    ReturnType<PluginUIContext['builders']['structure']['hierarchy']['applyPreset']>
>;

// Molecule interface representing a loaded molecule in the viewer.
export interface Molecule {
    label: string;
    name: string;
    filename: string;
    presetResult: PresetResult;
    trajectory: any;
    alignmentData?: AlignmentData | undefined;
}

/**
 * Loads a molecule file into the viewer.
 * @param plugin The PluginUIContext instance.
 * @param file The molecule file to load.
 * @param doGetAlignmentData Whether to extract alignment data.
 * @param centralise Whether to centralise coordinates.
 * @param alignment Optional alignment data to use.
 * @returns A promise that resolves to a Molecule object or undefined if loading fails.
 */
export async function loadMoleculeFileToViewer(
    plugin: PluginUIContext,
    file: Asset.File,
    doGetAlignmentData: boolean,
    centralise: boolean,
    alignment?: AlignmentData
): Promise<Molecule | undefined> {
    console.log('Loading molecule file into viewer:', file.name);
    const data = await plugin.builders.data.readFile(
        { file, label: file.name },
        { state: { isGhost: true } }
    );
    if (!data) {
        console.error('Failed to read file:', file.name);
        return;
    }
    console.log('File read successfully:', data);
    const myProvider = {
        ...RibocodeMmcifProvider,
        parse: async (
            plugin: PluginContext,
            data: StateObjectRef<PluginStateObject.Data.String | PluginStateObject.Data.Binary>,
            params: RibocodeMmcifParseParams | undefined
        ) => RibocodeMmcifProvider.parse(plugin, data,
            { ...params, centraliseCoordinates: centralise, alignmentData: alignment })
    };
    const trajectory = await plugin.builders.structure.parseTrajectory(data.data, myProvider);
    if (!trajectory) {
        console.error('Failed to parse trajectory from file:', file.name);
        return;
    }
    const presetResult = await plugin.builders.structure.hierarchy.applyPreset(
        trajectory, 'default');
    if (!presetResult) {
        console.error('Failed to apply preset to trajectory from file:', file.name);
        return;
    }
    console.log('applyPreset result:', presetResult);
    // Check presetResult and its representations
    // Check model.
    const model = presetResult.model;
    if (!model) {
        console.warn('No model found in presetResult.');
    }
    console.log('Model from presetResult:', model);
    // Check modelproperties.
    const modelProperties = presetResult.modelProperties;
    if (!modelProperties) {
        console.warn('No modelProperties found in presetResult.');
    }
    console.log('Model properties from presetResult:', modelProperties);
    // Check unitcell.
    const unitcell = presetResult.unitcell;
    if (!unitcell) {
        console.warn('No unitcell found in presetResult.');
    }
    console.log('Unitcell from presetResult:', unitcell);
    // Check representation structure.
    const structure = presetResult.structure;
    if (!structure) {
        console.warn('No structure found in presetResult.');
    }
    console.log('Structure from presetResult:', structure);
    // Check representation structureProperties.
    const structureProperties = presetResult.structureProperties;
    if (!structureProperties) {
        console.warn('No structureProperties found in presetResult.');
    }
    console.log('Structure properties from presetResult:', structureProperties);
    // Check representation.
    const representation = presetResult.representation;
    if (!representation) {
        console.warn('No representation found in presetResult.');
    }
    console.log('Representation from presetResult:', representation);
    // Check representation components.
    const components = representation.components;
    if (!components) {
        console.warn('No components found in representation.');
    }
    console.log('Components from representation:', components);
    // Check representation representations.
    const representations = representation.representations;
    if (!representations) {
        console.warn('No representations found in representation.');
    }
    console.log('Representations from representation:', representations);
    // Check for polymer representation
    const polymer = representations.polymer;
    if (!polymer) {
        console.warn('polymer representation was not created.');
    } else {
        console.log('polymer representation created:', polymer);
    }
    // Get alignment data if requested.
    let alignmentData: AlignmentData | undefined = undefined;
    if (doGetAlignmentData) {
        alignmentData = await getAlignmentData(plugin, trajectory);
    }
    const name: string = presetResult?.model.data?.label.toLocaleUpperCase() || file.name;
    const label: string = name.split(' ')[0];
    console.log('Molecule loaded:', label, name, file.name);
    return { label: label, name: name, filename: file.name, presetResult: presetResult, trajectory, alignmentData: alignmentData };
}

/**
 * Gets alignment data for a molecule including the type and location
 * of all atoms in the structure.
 * @param plugin The PluginUIContext instance.
 * @param trajectory The trajectory StateObjectRef.
 * @returns A promise that resolves to an object containing AlignmentData.
 */
export async function getAlignmentData(plugin: PluginUIContext, trajectory: any):
    Promise<AlignmentData> {
    console.log('Extracting alignment data');
    let x : number[] = [];
    let y : number[] = [];
    let z : number[] = [];
    let type : string[] = [];
    const trajCell = plugin.state.data.cells.get(trajectory.ref);
    if (trajCell?.obj?.data) {
        const frame = trajCell.obj.data.frames[0];
        console.log('Frame data:', frame);
        if (frame?.atomicConformation?.x) {
            x = Array.from(frame.atomicConformation.x);
        }
        if (frame?.atomicConformation?.y) {
            y = Array.from(frame.atomicConformation.y);
        }
        if (frame?.atomicConformation?.z) {
            z = Array.from(frame.atomicConformation.z);
        }
        if (frame?.atomicHierarchy?.atoms?.type_symbol?.__array) {
            type = frame.atomicHierarchy.atoms.type_symbol.__array as string[];
        }
    }
    // console.log('Alignment data extracted:', id.length, x.length, y.length, z.length, type.length);
    console.log('Alignment data extracted:', x.length, y.length, z.length, type.length);
    return { x, y, z, atomType: type };
}