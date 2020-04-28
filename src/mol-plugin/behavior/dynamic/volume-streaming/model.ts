/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { Volume } from '../../../../mol-model/volume';
import { Structure } from '../../../../mol-model/structure';

export class VolumeServerInfo extends PluginStateObject.Create<VolumeServerInfo.Data>({ name: 'Volume Streaming', typeClass: 'Object' }) { }

export namespace VolumeServerInfo {
    export type Kind = 'x-ray' | 'em'
    export interface EntryData  {
        kind: Kind,
        // for em, the EMDB access code, for x-ray, the PDB id
        dataId: string,
        header: VolumeServerHeader,
        emDefaultContourLevel?: Volume.IsoValue,
    }
    export interface Data {
        serverUrl: string,
        entries: EntryData[],
        structure: Structure
    }
}

export interface VolumeServerHeader {
    /** Format version number  */
    formatVersion: string,

    /** Axis order from the slowest to fastest moving, same as in CCP4 */
    axisOrder: number[],

    /** Origin in fractional coordinates, in axisOrder */
    origin: number[],

    /** Dimensions in fractional coordinates, in axisOrder */
    dimensions: number[],

    spacegroup: VolumeServerHeader.Spacegroup,
    channels: string[],

    /** Determines the data type of the values */
    valueType: VolumeServerHeader.ValueType,

    /** The value are stored in blockSize^3 cubes */
    blockSize: number,
    sampling: VolumeServerHeader.Sampling[],

    /** Precision data the server can show. */
    availablePrecisions: VolumeServerHeader.DetailLevel[],

    isAvailable: boolean
}

export namespace VolumeServerHeader {
    export type ValueType = 'float32' | 'int8'

    export namespace ValueType {
        export const Float32: ValueType = 'float32';
        export const Int8: ValueType = 'int8';
    }

    export type ValueArray = Float32Array | Int8Array

    export type DetailLevel = { precision: number, maxVoxels: number }

    export interface Spacegroup {
        number: number,
        size: number[],
        angles: number[],
        /** Determine if the data should be treated as periodic or not. (e.g. X-ray = periodic, EM = not periodic) */
        isPeriodic: boolean,
    }

    export interface ValuesInfo {
        mean: number,
        sigma: number,
        min: number,
        max: number
    }

    export interface Sampling {
        byteOffset: number,

        /** How many values along each axis were collapsed into 1 */
        rate: number,
        valuesInfo: ValuesInfo[],

        /** Number of samples along each axis, in axisOrder  */
        sampleCount: number[]
    }
}