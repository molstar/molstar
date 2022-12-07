/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ke Ma <mark.ma@rcsb.org>
 */
import { Structure } from '../../../mol-model/structure';
import { Vec3 } from '../../..//mol-math/linear-algebra/3d/vec3';
import { PluginContext } from '../../../mol-plugin/context';
import { CameraFocusOptions } from '../camera';
import { PrincipalAxes } from '../../../mol-math/linear-algebra/matrix/principal-axes';
import { StructureComponentRef } from '../structure/hierarchy-state';

export function getPolymerPositions(polymerStructure: Structure): Float32Array {
    let positionIndex = 0;
    const tmpMatrix = Vec3.zero();
    const positions = new Float32Array(polymerStructure.atomicResidueCount * 3);
    for (let i = 0; i < polymerStructure.units.length; i++) {
        const unit = polymerStructure.units[i];
        const { polymerElements } = unit.props;
        const readPosition = unit.conformation.position;
        if (polymerElements) {
            for (let j = 0; j < polymerElements.length; j++) {
                readPosition(polymerElements[j], tmpMatrix);
                positions[positionIndex] = tmpMatrix[0];
                positions[positionIndex + 1] = tmpMatrix[1];
                positions[positionIndex + 2] = tmpMatrix[2];
                positionIndex += 3;
            }
        }
    }
    return positions;
}
export function calculateDisplacement(position: Vec3, origin: Vec3, normalDir: Vec3) {
    const A = normalDir[0];
    const B = normalDir[1];
    const C = normalDir[2];
    const D = -A * origin[0] - B * origin[1] - C * origin[2];

    const x = position[0];
    const y = position[1];
    const z = position[2];

    const displacement = (A * x + B * y + C * z + D) / Math.sqrt(A * A + B * B + C * C);
    return displacement;
}

export function getAxesToFlip(position: Vec3, origin: Vec3, up: Vec3, normalDir: Vec3) {
    const toYAxis = calculateDisplacement(position, origin, normalDir);
    const toXAxis = calculateDisplacement(position, origin, up);
    const axes: ('aroundX' | 'aroundY')[] = [];
    if (toYAxis < 0) axes.push('aroundY');
    if (toXAxis < 0) axes.push('aroundX');
    return axes;
}

export function getFirstResidueOrAveragePosition(structure: Structure, caPositions: Float32Array): Vec3 {
    if (structure.units.length === 1) {
        // if only one chain => first residue coordinates
        return Vec3.create(caPositions[0], caPositions[1], caPositions[2]);
    } else {
        // if more than one chain => average of coordinates of the first chain
        const tmpMatrixPos = Vec3.zero();
        const atomIndices = structure.units[0].props.polymerElements;
        const firstChainPositions = [];
        if (atomIndices) {
            for (let i = 0; i < atomIndices.length; i++) {
                const coordinates = structure.units[0].conformation.position(atomIndices[i], tmpMatrixPos);
                for (let j = 0; j < coordinates.length; j++) {
                    firstChainPositions.push(coordinates[j]);
                }
            }
            let sumX = 0;
            let sumY = 0;
            let sumZ = 0;
            for (let i = 0; i < firstChainPositions.length; i += 3) {
                sumX += firstChainPositions[i];
                sumY += firstChainPositions[i + 1];
                sumZ += firstChainPositions[i + 2];
            }
            const averagePosition = Vec3.zero();
            averagePosition[0] = sumX / atomIndices.length;
            averagePosition[1] = sumY / atomIndices.length;
            averagePosition[2] = sumZ / atomIndices.length;
            return averagePosition;
        } else {
            return Vec3.create(caPositions[0], caPositions[1], caPositions[2]);
        }
    }

}

export function pcaFocus(plugin: PluginContext, options: Partial<CameraFocusOptions> & { principalAxes?: PrincipalAxes, positionToFlip?: Vec3 }) {
    if (options?.principalAxes) {
        const { origin, dirA, dirB, dirC } = options.principalAxes.boxAxes;
        let toFlip: ('aroundX' | 'aroundY')[] = [];
        if (options.positionToFlip) {
            toFlip = getAxesToFlip(options.positionToFlip, origin, dirA, dirB);
        }
        toFlip.forEach((axis)=>{
            if (axis === 'aroundY') {
                Vec3.negate(dirC, dirC);
            } else if (axis === 'aroundX') {
                Vec3.negate(dirA, dirA);
                Vec3.negate(dirC, dirC);
            }
        });
        if (plugin.canvas3d) {
            const position = Vec3();
            Vec3.scaleAndAdd(position, position, origin, 100);
            plugin.canvas3d.camera.setState({ position }, 0);
            const deltaDistance = Vec3();
            Vec3.negate(deltaDistance, position);
            if (Vec3.dot(deltaDistance, dirC) <= 0) {
                Vec3.negate(plugin.canvas3d.camera.position, position);
            }
            const up = Vec3.create(0, 1, 0);
            if (Vec3.dot(up, dirA) <= 0) {
                Vec3.negate(plugin.canvas3d?.camera.up, plugin.canvas3d.camera.up);
            }
        }
        return { origin, dirA, dirB, dirC };
    }
    return {
        origin: Vec3.zero(),
        dirA: Vec3.zero(),
        dirB: Vec3.zero(),
        dirC: Vec3.zero()
    };
}

export function getPcaTransform(group: StructureComponentRef[]): { principalAxes?: PrincipalAxes, positionToFlip?: Vec3 } | undefined {
    const polymerStructure = group[0].cell.obj?.data;
    if (!polymerStructure) {
        return undefined;
    }
    // if ('_pcaTransformData' in polymerStructure.currentPropertyData) {
    //     console.log("run the cache")
    //     console.log(polymerStructure.currentPropertyData._pcaTransformData)
    //     return polymerStructure.currentPropertyData._pcaTransformData;
    // }
    if (!polymerStructure.units[0]?.props.polymerElements?.length) {
        polymerStructure.currentPropertyData._pcaTransformData = undefined;
        return undefined;
    }
    const positions = getPolymerPositions(polymerStructure);
    const positionToFlip = getFirstResidueOrAveragePosition(polymerStructure, positions);
    polymerStructure.currentPropertyData._pcaTransformData = { principalAxes: PrincipalAxes.ofPositions(positions), positionToFlip };
    return { principalAxes: PrincipalAxes.ofPositions(positions), positionToFlip };
}
