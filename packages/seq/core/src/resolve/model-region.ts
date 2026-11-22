/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isColumn, type Node, type Value } from '@molstar/core/repr/tree';
import type { Coverage, MatrixItems } from '@molstar/seqx-core/data';
import { channelAddress } from '@molstar/seqx-core/data/channels';
import { inlinePlan, sourcePlan } from './model-entries.js';
import type { Placed } from './model-molecules.js';
import type { AnnotationEntry, SystemEntry } from './model.js';

// Pairs, grids and matrices read by region (_spec/README.md "Resolve", Access plans by source): a `pairs` or `grid`
// node with a `source` reads it by region through the store on its first system (`pairs.a`, `grid.x`), its second
// system the entry's `ySystem`; its `field` a channel of the source. A `matrix` node is on its `target`, its states a
// system of their own (`<matrix>/(states)`, a position a state, EC-35's rows), their count from the inline `states` or,
// once the source is open, its coverage (`ShapeInfo.states`, as a length: `sourceLengths`).

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

/** A matrix's states system's address: its node's, then `(states)`. */
export function statesSystemOf(matrix: string): string {
  return `${matrix}/(states)`;
}

/**
 * The member system an address is (an alignment's rows, `<alignment>/(rows)`; a matrix's states, `<matrix>/(states)`:
 * a position a named member, EC-35): its owner node's address and which; undefined for any other.
 */
export function memberSystem(system: string): { owner: string; kind: 'rows' | 'states' } | undefined {
  if (system.endsWith('/(rows)')) return { owner: system.slice(0, -7), kind: 'rows' };
  if (system.endsWith('/(states)')) return { owner: system.slice(0, -9), kind: 'states' };
  return undefined;
}

/** The name a matrix's states system has (what no source names: its states are its coverage's `states`). */
export const StatesSystemName = '(states)';

/** The source a node reads: its channel when it names a `field`. */
function readOf(node: Node, source: string): string {
  const field = p(node, 'field');
  return typeof field === 'string' ? channelAddress(source, field) : source;
}

/** A `pairs` or `grid` node with a `source`: read by region on its first system, its second the entry's `ySystem`. */
export function regionEntry(at: Placed, x: SystemEntry, y: SystemEntry, source: string): AnnotationEntry {
  const shape = at.node.kind as 'pairs' | 'grid';
  return {
    address: at.address,
    shape,
    target: x.address,
    ySystem: y.address,
    fields: shape === 'grid' ? ['value'] : [],
    access: sourcePlan(readOf(at.node, source), x.name, shape),
  };
}

/**
 * A `matrix` node: its states system (its length the inline `states`' count, else `knownStates` from the source's
 * coverage, else unknown: pending) and its annotation, inline cells an items table of one item a position of `target`.
 */
export function matrixEntry(
  at: Placed,
  target: SystemEntry,
  source: string | undefined,
  knownStates: number | undefined,
): { states: SystemEntry; annotation: AnnotationEntry } {
  const names = p(at.node, 'states');
  const inlineStates = Array.isArray(names) ? (names as string[]) : undefined;
  const count = inlineStates?.length ?? knownStates;
  const statesAddress = statesSystemOf(at.address);
  const states: SystemEntry = {
    address: statesAddress,
    kind: 'sequence',
    ...(count !== undefined ? { length: count } : {}),
    cyclic: false,
    name: StatesSystemName,
    numberings: ['index', 'ordinal'],
    // the states are its source's coverage's, never letters to read
    access: sourcePlan(source ?? '', StatesSystemName, 'matrix'),
  };
  const base = {
    address: at.address,
    shape: 'matrix' as const,
    target: target.address,
    fields: [],
    states: statesAddress,
  };
  if (source !== undefined)
    return { states, annotation: { ...base, access: sourcePlan(readOf(at.node, source), target.name, 'matrix') } };
  const values = p(at.node, 'values');
  const codes = p(at.node, 'codes');
  const cells = isColumn(values) ? (values as Float32Array) : isColumn(codes) ? (codes as Uint8Array) : undefined;
  const n = cells && count ? Math.floor(cells.length / count) : 0;
  const matrix: MatrixItems = {
    form: 'items',
    start: Int32Array.from({ length: n }, (_, i) => i),
    end: Int32Array.from({ length: n }, (_, i) => i + 1),
    states: count ?? 0,
    ...(cells instanceof Float32Array ? { value: cells.slice(0, n * (count ?? 0)) } : {}),
    ...(cells instanceof Uint8Array ? { codes: cells.slice(0, n * (count ?? 0)) } : {}),
  };
  const length = target.length ?? n;
  return {
    states,
    annotation: { ...base, access: inlinePlan(at.address, { ...target, length }, 'matrix', { matrix }) },
  };
}

/** The states a matrix source's coverage lists for its target system (by name), when it lists them. */
export function coverageStates(coverage: Coverage | undefined, target: SystemEntry): readonly string[] | undefined {
  return coverage?.shapes.find((s) => s.shape === 'matrix' && s.system === target.name)?.states;
}
