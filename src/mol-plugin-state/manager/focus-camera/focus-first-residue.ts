/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ke Ma <mark.ma@rcsb.org>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
import { Structure } from '../../../mol-model/structure';
import { Vec3 } from '../../..//mol-math/linear-algebra/3d/vec3';
import { PluginContext } from '../../../mol-plugin/context';
import { PrincipalAxes } from '../../../mol-math/linear-algebra/matrix/principal-axes';
import { StructureComponentRef } from '../structure/hierarchy-state';
import { Camera } from '../../../mol-canvas3d/camera';

function getPolymerResiduePositions(structure: Structure): Float32Array | undefined {
    if (structure.atomicResidueCount === 1) return undefined;

    const { polymerResidueCount } = structure;
    if (polymerResidueCount <= 1) return undefined;

    const stride = 2 ** Math.max(Math.ceil(Math.log10(polymerResidueCount / 1000)), 0);
    const size = stride === 1
        ? polymerResidueCount
        : Math.ceil(polymerResidueCount / stride) + structure.units.length;

    const tmpPos = Vec3();
    const positions = new Float32Array(3 * size);
    let o = 0;
    for (const unit of structure.units) {
        const { polymerElements } = unit.props;
        const { conformation } = unit;
        if (polymerElements) {
            for (let i = 0; i < polymerElements.length; i += stride) {
                conformation.position(polymerElements[i], tmpPos);
                Vec3.toArray(tmpPos, positions, 3 * o);
                o++;
            }
        }
    }

    return positions.length !== o ? positions.slice(0, 3 * o) : positions;
}

function calculateDisplacement(position: Vec3, origin: Vec3, normalDir: Vec3) {
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

function getAxesToFlip(position: Vec3, origin: Vec3, up: Vec3, normalDir: Vec3) {
    const toYAxis = calculateDisplacement(position, origin, normalDir);
    const toXAxis = calculateDisplacement(position, origin, up);
    return {
        aroundX: toXAxis < 0,
        aroundY: toYAxis < 0,
    };
}

function getFirstResidueOrAveragePosition(structure: Structure, polymerPositions: Float32Array): Vec3 {
    if (structure.units.length === 1) {
        // if only one chain, return the coordinates of the first residue
        return Vec3.create(polymerPositions[0], polymerPositions[1], polymerPositions[2]);
    } else {
        // if more than one chain, return average of the coordinates of the first polymer chain
        const firstPolymerUnit = structure.units.find(u => u.props.polymerElements);
        if (firstPolymerUnit) {
            const pos = Vec3();
            const center = Vec3();
            const { polymerElements, conformation } = firstPolymerUnit;
            for (let i = 0, il = polymerElements.length; i < il; i++) {
                conformation.position(polymerElements[i], pos);
                Vec3.add(center, center, pos);
            }
            return Vec3.scale(center, center, 1 / polymerElements.length);
        } else {
            return Vec3.create(polymerPositions[0], polymerPositions[1], polymerPositions[2]);
        }
    }
}

export function pcaFocus(plugin: PluginContext, radius: number, options: { principalAxes: PrincipalAxes, positionToFlip?: Vec3 }) {
    if (!plugin.canvas3d) return;

    const { origin, dirA, dirB, dirC } = options.principalAxes.boxAxes;

    const up = Vec3.clone(dirA);
    const dir = Vec3.clone(dirC);

    if (options.positionToFlip) {
        const { aroundX, aroundY } = getAxesToFlip(options.positionToFlip, origin, up, dirB);

        if (aroundX) {
            Vec3.negate(dir, dir);
            Vec3.negate(up, up);
        }
        if (aroundY) {
            Vec3.negate(dir, dir);
        }
    }

    const position = Vec3.scale(Vec3(), origin, -100);
    if (Vec3.dot(position, up) <= 0) {
        Vec3.negate(dir, dir);
    }
    if (Vec3.dot(Vec3.unitY, dir) <= 0) {
        Vec3.negate(up, up);
    }

    return plugin.canvas3d.camera.getFocus(origin, radius, up, dir, Camera.createDefaultSnapshot());
}

interface PcaTransformData {
    principalAxes: PrincipalAxes;
    positionToFlip: Vec3;
}

export function getPcaTransform(group: StructureComponentRef[]): PcaTransformData | undefined {
    const structure = group[0].cell.obj?.data;
    if (!structure) return undefined;

    if ('_pcaTransformData' in structure.currentPropertyData) {
        return structure.currentPropertyData._pcaTransformData;
    }
    const positions = getPolymerResiduePositions(structure);
    if (!positions) {
        structure.currentPropertyData._pcaTransformData = undefined;
        return undefined;
    }
    const positionToFlip = getFirstResidueOrAveragePosition(structure, positions);
    const pcaTransfromData: PcaTransformData = {
        principalAxes: PrincipalAxes.ofPositions(positions),
        positionToFlip
    };
    structure.currentPropertyData._pcaTransformData = pcaTransfromData;
    return pcaTransfromData;
}
