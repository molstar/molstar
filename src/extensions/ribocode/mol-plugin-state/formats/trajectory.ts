import { StateTransforms } from '../../../../mol-plugin-state/transforms';
import { TrajectoryFormatCategory, TrajectoryFormatProvider } from '../../../../mol-plugin-state/formats/trajectory';
import { guessCifVariant } from '../../../../mol-plugin-state/formats/provider';
import { StateObjectRef } from '../../../../mol-state';
import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { PluginContext } from '../../../../mol-plugin/context';
import { computeBigMean, alignDataset } from '../../utils/geometry';

/**
 * Gets the default visuals for the trajectory.
 * @param plugin The PluginContext.
 * @param data The trajectory data.
 * @returns The default visuals for the trajectory.
 */
function defaultVisuals(plugin: PluginContext, data: { trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory> }) {
    return plugin.builders.structure.hierarchy.applyPreset(data.trajectory, 'default');
}

/**
 * Parameters for parsing mmCIF files with Ribocode extensions.
 * @interface RibocodeMmcifParseParams
 * @property {string | string[]} [trajectoryTags] - Tags to apply to the trajectory.
 * @property {boolean} [centraliseCoordinates] - Whether to centralise coordinates.
 * @property {any} [alignmentData] - Data for aligning the coordinates.
 */
export interface RibocodeMmcifParseParams {
    trajectoryTags?: string | string[];
    centraliseCoordinates?: boolean;
    alignmentData?: any;
}

/**
 * Ribocode mmCIF format provider for parsing mmCIF files with Ribocode-specific extensions.
 * @constant RibocodeMmcifProvider
 * @type {TrajectoryFormatProvider<RibocodeMmcifParseParams>}
 * @description A TrajectoryFormatProvider for mmCIF files with Ribocode extensions.
 */
export const RibocodeMmcifProvider: TrajectoryFormatProvider<RibocodeMmcifParseParams> = {
    label: 'mmCIF',
    description: 'mmCIF',
    category: TrajectoryFormatCategory,
    stringExtensions: ['cif', 'mmcif', 'mcif'],
    binaryExtensions: ['bcif'],
    isApplicable: (info, data) => {
        if (info.ext === 'mmcif' || info.ext === 'mcif') return true;
        // assume undetermined cif/bcif files are mmCIF
        if (info.ext === 'cif' || info.ext === 'bcif') return guessCifVariant(info, data) === -1;
        return false;
    },
    parse: async (plugin, data, params?: RibocodeMmcifParseParams) => {
        console.log('MmcifProvider.parse called');
        const state = plugin.state.data;
        console.log('Building state tree for data:', data);
        const cif = state.build().to(data)
            .apply(StateTransforms.Data.ParseCif, void 0, { state: { isGhost: true } });
        console.log('cif state built:', cif);
        console.log('Applying TrajectoryFromMmCif transform');
        const trajectory = await cif
            .apply(StateTransforms.Model.TrajectoryFromMmCif, void 0, { tags: params?.trajectoryTags })
            .commit({ revertOnError: true });
        if (params?.centraliseCoordinates || params?.alignmentData) {
            // Step 1: Get the cell from the state using trajectory.ref
            const trajCell = plugin.state.data.cells.get(trajectory.ref);
            // Step 2: Inspect the cell and its data
            // console.log('trajCell:', trajCell);
            // console.log('trajCell.obj:', trajCell?.obj);
            if (trajCell?.obj?.data) {
                // console.log('trajCell.obj.data:', trajCell?.obj?.data);
                // Inspect the first frame and representative
                const nframes = trajCell.obj.data.frames.length;
                if (nframes === 0) {
                    console.warn('No frames found in trajectory data.');
                    return { trajectory };
                } else if (nframes > 1) {
                    console.warn(`Multiple frames (${nframes}) found in trajectory data. Centralisation/alignment will be applied only to the first frame.`);
                }
                // console.log('Representative:', trajCell.obj.data.representative);
                const frame = trajCell.obj.data.frames[0];
                if (frame) {
                    // console.log('Frame data:', frame);
                    if (frame.atomicConformation && frame.atomicHierarchy) {
                        // console.log('atomicConformation:', frame.atomicConformation);
                        const x: number[] = Array.from(frame.atomicConformation.x);
                        const y: number[] = Array.from(frame.atomicConformation.y);
                        const z: number[] = Array.from(frame.atomicConformation.z);
                        // console.log('atomicHierarchy:', frame.atomicHierarchy);
                        const type_symbol = frame.atomicHierarchy.atoms.type_symbol.__array;
                        // console.log('Atom coordinates and type:', { x, y, z, type_symbol });

                        // Count of all different atom types
                        const atomTypeCounts: { [key: string]: number } = {};
                        for (const type of type_symbol) {
                            atomTypeCounts[type] = (atomTypeCounts[type] || 0) + 1;
                        }
                        console.log('Atom type counts:', atomTypeCounts);

                        const n = x.length;
                        const d = Math.floor(n / 5);
                        let newX: number[] = [];
                        let newY: number[] = [];
                        let newZ: number[] = [];
                        console.log('Centralising coordinates');
                        // Centralise coordinates logic
                        // 1. Compute centroid
                        const xMean: number = computeBigMean(x);
                        const yMean: number = computeBigMean(y);
                        const zMean: number = computeBigMean(z);
                        console.log('xMean:', xMean, 'yMean:', yMean, 'zMean:', zMean);
                        // 3. Calculate new coordinates
                        // Subtract centroid from each coordinate
                        for (let i = 0; i < n; i++) {
                            newX.push(x[i] - xMean);
                            newY.push(y[i] - yMean);
                            newZ.push(z[i] - zMean);
                            if (i % d === 0) {
                                console.log(`Centralised coordinates for atom ${i} ${type_symbol[i]} from (${x[i]}, ${y[i]}, ${z[i]}) to (${newX[i]}, ${newY[i]}, ${newZ[i]})`);
                            }
                        }
                        if (params?.alignmentData) {
                            // Use alignmentData to align the new coordinates
                            // console.log('Alignment data provided:', params.alignmentData);
                            const alignedCoordinates = alignDataset( type_symbol,
                                newX, newY, newZ,
                                params.alignmentData.type, params.alignmentData.x, params.alignmentData.y, params.alignmentData.z);
                            // Update newX, newY, newZ with aligned coordinates
                            newX = alignedCoordinates.alignedX;
                            newY = alignedCoordinates.alignedY;
                            newZ = alignedCoordinates.alignedZ;
                            console.log('Coordinates aligned using provided alignment data. ', newX.length);
                        } else {
                            console.log('No alignment data provided.');
                        }
                        // Replace coordinates
                        frame.atomicConformation.x = newX;
                        frame.atomicConformation.y = newY;
                        frame.atomicConformation.z = newZ;
                    } else {
                        // Log all keys to help discover coordinate storage
                        console.log('Frame keys:', Object.keys(frame));
                    }
                } else {
                    console.log('No frame data found.');
                }
            } else {
                console.log('trajCell or trajCell.obj.data is undefined.');
            }
        } else {
            console.log('Centralise coordinates not requested.');
        }
        console.log('TrajectoryFromMmCif applied, trajectory:', trajectory);
        console.log('MmcifProvider.parse completed');
        if ((cif.selector.cell?.obj?.data.blocks.length || 0) > 1) {
            plugin.state.data.updateCellState(cif.ref, { isGhost: false });
        }
        return { trajectory };
    },
    visuals: defaultVisuals
};