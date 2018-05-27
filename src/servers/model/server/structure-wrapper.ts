/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from 'mol-model/structure';

export enum StructureSourceType {
    File,
    Cache
}

export interface StructureInfo {
    sourceType: StructureSourceType;
    readTime: number;
    parseTime: number;

    sourceId: string,
    entryId: string,
    filename: string
}

export class StructureWrapper {
    info: StructureInfo;

    key: string;
    approximateSize: number;
    structure: Structure;
}

export function getStructure(filename: string): Promise<StructureWrapper>
export function getStructure(sourceId: string, entryId: string): Promise<StructureWrapper>
export function getStructure(sourceIdOrFilename: string, entryId?: string): Promise<StructureWrapper> {
    return 0 as any;
}