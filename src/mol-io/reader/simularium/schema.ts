/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Simularium trajectory file format (Allen Institute for Cell Science).
 * https://github.com/simularium/simulariumio
 */

/** A single entry of the `typeMapping` describing one agent type. */
export interface SimulariumTypeEntry {
    readonly name: string
    readonly geometry?: {
        /** e.g. `SPHERE`, `FIBER`, `PDB`, `OBJ`, `SPHERE_GROUP`. */
        readonly displayType?: string
        readonly url?: string
        readonly color?: string
    }
}

/** Physical unit (`magnitude` of `name`, e.g. 0.1 of `nm`). */
export interface SimulariumUnit {
    readonly magnitude: number
    readonly name: string
}

export interface SimulariumTrajectoryInfo {
    readonly version: number
    readonly totalSteps: number
    readonly timeStepSize?: number
    readonly spatialUnits?: SimulariumUnit
    readonly timeUnits?: SimulariumUnit
    readonly size?: { readonly x: number, readonly y: number, readonly z: number }
    readonly trajectoryTitle?: string
    /** Map from agent type id (as a string key) to its metadata. */
    readonly typeMapping: { readonly [id: string]: SimulariumTypeEntry }
}

export interface SimulariumFrame {
    readonly frameNumber: number
    readonly time: number
    /**
     * Flat agent buffer for this frame, packed as a sequence of agents. Each agent is
     * `SimulariumMinValuesPerAgent` fixed floats followed by `nSubpoints` subpoint floats
     * (see `SimulariumAgentBuffer`). The buffer starts at the first agent's `visType`,
     * i.e. it does NOT include the per-frame `frameNumber`/`time`/`nAgents` prefix.
     */
    readonly data: Float32Array
}

export interface SimulariumFile {
    readonly trajectoryInfo: SimulariumTrajectoryInfo
    readonly frames: ReadonlyArray<SimulariumFrame>
}

/** Float offsets of the fixed fields within a single agent in `SimulariumFrame.data`. */
export const SimulariumAgentBuffer = {
    VIS_TYPE: 0,
    INSTANCE_ID: 1,
    TYPE_ID: 2,
    POS_X: 3,
    POS_Y: 4,
    POS_Z: 5,
    ROT_X: 6,
    ROT_Y: 7,
    ROT_Z: 8,
    /** Collision radius. */
    RADIUS: 9,
    /** Number of subpoint floats that follow the fixed fields. */
    N_SUBPOINTS: 10,
} as const;

/** Number of fixed floats per agent before the variable-length subpoints. */
export const SimulariumMinValuesPerAgent = 11;

/** Values of the agent `visType` field. */
export const SimulariumVisType = {
    DEFAULT: 1000,
    FIBER: 1001,
} as const;
