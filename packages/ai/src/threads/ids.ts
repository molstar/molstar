/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

declare const threadIdBrand: unique symbol;
declare const threadEntryIdBrand: unique symbol;
declare const threadOwnerIdBrand: unique symbol;
declare const threadCommandIdBrand: unique symbol;

export type ThreadId = string & { readonly [threadIdBrand]: undefined };
export type ThreadEntryId = string & { readonly [threadEntryIdBrand]: undefined };
export type ThreadOwnerId = string & { readonly [threadOwnerIdBrand]: undefined };
export type ThreadCommandId = string & { readonly [threadCommandIdBrand]: undefined };

function nonempty(value: string, label: string): string {
  if (value.trim().length === 0) throw new Error(`${label} must be nonempty`);
  return value;
}

export function threadId(value: string): ThreadId {
  return nonempty(value, 'ThreadId') as ThreadId;
}

export function threadEntryId(value: string): ThreadEntryId {
  return nonempty(value, 'ThreadEntryId') as ThreadEntryId;
}

export function threadOwnerId(value: string): ThreadOwnerId {
  return nonempty(value, 'ThreadOwnerId') as ThreadOwnerId;
}

export function threadCommandId(value: string): ThreadCommandId {
  return nonempty(value, 'ThreadCommandId') as ThreadCommandId;
}
