/**
 * Copyright (c) 2022-23 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ke Ma <mark.ma@rcsb.org>
 * @author David Sehnal <david.sehnal@gmail.com>
 */
import { Structure } from '../../../mol-model/structure';
import { Vec3 } from '../../..//mol-math/linear-algebra/3d/vec3';
import { PluginContext } from '../../../mol-plugin/context';
import { PrincipalAxes } from '../../../mol-math/linear-algebra/matrix/principal-axes';
import { StructureComponentRef } from '../structure/hierarchy-state';
import { Camera } from '../../../mol-canvas3d/camera';


function getPolymerPositions(structure: Structure): Float32Array | undefined {
    if (structure.atomicResidueCount === 1) return undefined;

    let polymerElementCount = 0;
    for (const unit of structure.units) {
        polymerElementCount += unit.props.polymerElements?.length ?? 0;
    }

    if (polymerElementCount <= 1) return undefined;

    const stride = 2 ** Math.max(Math.ceil(Math.log10(polymerElementCount / 1000)), 0);
    const size = stride === 1
        ? polymerElementCount
        : Math.ceil(polymerElementCount / stride) + structure.units.length;

    const tmpPos = Vec3.zero();
    const positions = new Float32Array(3 * size);
    let o = 0;
    for (const unit of structure.units) {
        const { polymerElements } = unit.props;
        const { position } = unit.conformation;
        if (polymerElements) {
            for (let i = 0; i < polymerElements.length; i += stride) {
                position(polymerElements[i], tmpPos);
                Vec3.toArray(tmpPos, positions, 3 * o);
                o++;
            }
        }
    }
    if (positions.length !== o) return positions.slice(0, 3 * o);
    return positions;
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
        // if only one chain => first residue coordinates
        return Vec3.create(polymerPositions[0], polymerPositions[1], polymerPositions[2]);
    } else {
        // if more than one chain => average of coordinates of the first chain
        const polymerElements = structure.units.find(u => u.props.polymerElements)?.props.polymerElements;
        if (polymerElements?.length) {
            const pos = Vec3.zero();
            const center = Vec3.zero();
            const { position } = structure.units[0].conformation;
            for (let i = 0; i < polymerElements.length; i++) {
                position(polymerElements[i], pos);
                Vec3.add(center, center, pos);
            }
            Vec3.scale(center, center, 1 / polymerElements.length);
            return center;
        } else {
            return Vec3.create(polymerPositions[0], polymerPositions[1], polymerPositions[2]);
        }
    }
}

export function pcaFocus(plugin: PluginContext, radius: number, options: { principalAxes: PrincipalAxes, positionToFlip?: Vec3 }) {
    const { origin, dirB } = options.principalAxes.boxAxes;
    let { dirA: up, dirC: dir } = options.principalAxes.boxAxes;

    if (options.positionToFlip) {
        const { aroundX, aroundY } = getAxesToFlip(options.positionToFlip, origin, up, dirB);

        // Clone the up and dir since we will be mutating them below
        up = Vec3.clone(up);
        dir = Vec3.clone(dir);

        if (aroundX) {
            Vec3.negate(dir, dir);
            Vec3.negate(up, up);
        }
        if (aroundY) {
            Vec3.negate(dir, dir);
        }
    }

    if (plugin.canvas3d) {
        const position = Vec3();
        // NOTE: the below Vec3.scale is simplification of
        //   Vec3.scaleAndAdd(position, position, origin, 100);
        //   plugin.canvas3d.camera.setState({ position }, 0);
        //   const deltaDistance = Vec3();
        //   Vec3.negate(deltaDistance, position);
        // from the original code.
        Vec3.scale(position, origin, -100);
        if (Vec3.dot(position, up) <= 0) {
            Vec3.negate(dir, dir);
        }
        const upY = Vec3.create(0, 1, 0);
        if (Vec3.dot(upY, dir) <= 0) {
            Vec3.negate(up, up);
        }
    }

    return plugin.canvas3d?.camera.getFocus(origin, radius, up, dir, Camera.createDefaultSnapshot());
}

export interface PCAFocusInfo {
    principalAxes: PrincipalAxes;
    positionToFlip: Vec3;
}

export function getPcaTransform(group: StructureComponentRef[]): PCAFocusInfo | undefined {
    const structure = group[0].cell.obj?.data;
    if (!structure) {
        return undefined;
    }
    if ('_pcaTransformData' in structure.currentPropertyData) {
        return structure.currentPropertyData._pcaTransformData;
    }
    const positions = getPolymerPositions(structure);
    if (!positions) {
        structure.currentPropertyData._pcaTransformData = undefined;
        return undefined;
    }
    const positionToFlip = getFirstResidueOrAveragePosition(structure, positions);
    structure.currentPropertyData._pcaTransformData = { principalAxes: PrincipalAxes.ofPositions(positions), positionToFlip } as PCAFocusInfo;
    return structure.currentPropertyData._pcaTransformData;
}
