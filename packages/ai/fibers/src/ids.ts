/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

declare const fiberIdBrand: unique symbol;
declare const fiberEntryIdBrand: unique symbol;
declare const fiberOwnerIdBrand: unique symbol;
declare const fiberCommandIdBrand: unique symbol;

export type FiberId = string & { readonly [fiberIdBrand]: undefined };
export type FiberEntryId = string & { readonly [fiberEntryIdBrand]: undefined };
export type FiberOwnerId = string & { readonly [fiberOwnerIdBrand]: undefined };
export type FiberCommandId = string & { readonly [fiberCommandIdBrand]: undefined };

function nonempty(value: string, label: string): string {
  if (value.trim().length === 0) throw new Error(`${label} must be nonempty`);
  return value;
}

export function fiberId(value: string): FiberId {
  return nonempty(value, 'FiberId') as FiberId;
}

export function fiberEntryId(value: string): FiberEntryId {
  return nonempty(value, 'FiberEntryId') as FiberEntryId;
}

export function fiberOwnerId(value: string): FiberOwnerId {
  return nonempty(value, 'FiberOwnerId') as FiberOwnerId;
}

export function fiberCommandId(value: string): FiberCommandId {
  return nonempty(value, 'FiberCommandId') as FiberCommandId;
}
