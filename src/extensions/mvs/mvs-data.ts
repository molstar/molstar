/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

// Re-export everything from the molviewspec package
export * from "molviewspec";

// Import types and values separately for proper re-export
import {
  MVSData as _MVSData,
  MVSTreeSchema as _MVSTreeSchema,
  createMVSBuilder as _createMVSBuilder,
  GlobalMetadata as _GlobalMetadata,
} from "molviewspec";

import type {
  MVSData as MVSDataType,
  MVSTree as MVSTreeType,
  SnapshotMetadata as SnapshotMetadataType,
  Snapshot as SnapshotType,
  MVSData_State as MVSData_StateType,
  MVSData_States as MVSData_StatesType,
  Root as RootType,
} from "molviewspec";

// Re-export values for compatibility
export const MVSData = _MVSData;
export const MVSTreeSchema = _MVSTreeSchema;
export const createMVSBuilder = _createMVSBuilder;
export const GlobalMetadata = _GlobalMetadata;

// Re-export types
export type MVSData = MVSDataType;
export type MVSTree = MVSTreeType;
export type SnapshotMetadata = SnapshotMetadataType;
export type Snapshot = SnapshotType;
export type MVSData_State = MVSData_StateType;
export type MVSData_States = MVSData_StatesType;
export type Root = RootType;
