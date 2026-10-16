/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { JsonObject, ModelVisiblePart, Reference } from '../../messages/index.js';
import { UnrepresentableInput } from '../../messages/index.js';

const strictUtf8 = new TextDecoder('utf-8', { fatal: true, ignoreBOM: true });

export type AttachmentPlacement = 'system' | 'user' | 'assistant' | 'tool-result';

/** Retained, JSON-compatible immutable inline bytes. Never a fetchable URL. */
export interface InlineAttachment {
  readonly mediaType: string;
  readonly data: string;
  readonly byteLength: number;
  readonly digest: string;
  readonly name?: string;
}

function mediaType(value: string): string {
  const normalized = value.trim().toLowerCase();
  if (!/^[a-z0-9!#$&^_.+-]+\/[a-z0-9!#$&^_.+-]+$/.test(normalized)) {
    throw new Error('Inline attachment media type must be an explicit MIME type.');
  }
  return normalized;
}

export function decodeBase64(data: string): Uint8Array {
  if (!/^(?:[A-Za-z0-9+/]{4})*(?:[A-Za-z0-9+/]{2}==|[A-Za-z0-9+/]{3}=)?$/.test(data)) {
    throw new Error('Inline attachment has malformed base64 encoding.');
  }
  const decoded = atob(data);
  if (btoa(decoded) !== data) throw new Error('Inline attachment has non-canonical base64 encoding.');
  return Uint8Array.from(decoded, (character) => character.charCodeAt(0));
}

function encodeBase64(bytes: Uint8Array): string {
  let binary = '';
  // Avoid argument-spreading large attachments onto the JS call stack.
  for (const byte of bytes) binary += String.fromCharCode(byte);
  return btoa(binary);
}

/** `sha256:<lowercase hex>` of inline bytes, the form a reference's `digest` takes. */
export async function bytesDigest(bytes: Uint8Array): Promise<string> {
  const hash = new Uint8Array(await crypto.subtle.digest('SHA-256', new Uint8Array(bytes)));
  return `sha256:${Array.from(hash, (byte) => byte.toString(16).padStart(2, '0')).join('')}`;
}

function decodePercent(body: string): Uint8Array {
  const decoded: number[] = [];
  for (let index = 0; index < body.length; index++) {
    const character = body[index] ?? '';
    if (character === '%') {
      const hex = body.slice(index + 1, index + 3);
      if (!/^[a-f0-9]{2}$/i.test(hex)) throw new Error('Inline attachment has malformed percent encoding.');
      decoded.push(Number.parseInt(hex, 16));
      index += 2;
    } else {
      const byte = character.charCodeAt(0);
      if (byte < 0x21 || byte > 0x7e || character === '#' || character === '?') {
        throw new Error('Inline attachment data URI must percent-encode non-URI bytes.');
      }
      decoded.push(byte);
    }
  }
  return Uint8Array.from(decoded);
}

/** Recognize all data URIs before consulting provider support, including unsupported MIME types.
 * Malformed or mismatched retained bytes are plain errors: the input is corrupt, not unrepresentable. */
export async function inlineAttachment(reference: Reference): Promise<InlineAttachment | undefined> {
  if (!/^data:/i.test(reference.uri)) return undefined;
  const comma = reference.uri.indexOf(',');
  if (comma < 0) throw new Error('Inline attachment data URI is missing its comma separator.');
  const fields = reference.uri.slice(5, comma).split(';');
  const embeddedType = mediaType(fields.shift() || 'text/plain');
  let base64 = false;
  let charset: string | undefined;
  for (const [index, field] of fields.entries()) {
    if (field.toLowerCase() === 'base64' && index === fields.length - 1) {
      base64 = true;
    } else if (!/^charset=(?:utf-8|us-ascii)$/i.test(field) || !embeddedType.startsWith('text/')) {
      throw new Error('Inline attachment data URI has unsupported encoding parameters.');
    } else {
      if (charset !== undefined) throw new Error('Inline attachment data URI repeats charset parameters.');
      charset = field.slice(8).toLowerCase();
    }
  }
  if (reference.mediaType !== undefined) {
    const declaredFields = reference.mediaType.split(';');
    if (mediaType(declaredFields[0] ?? '') !== embeddedType) {
      throw new Error('Inline attachment declared media type does not match its data URI.');
    }
    if (
      (declaredFields.length > 1 && !embeddedType.startsWith('text/')) ||
      declaredFields.length > 2 ||
      declaredFields.slice(1).some((field) => !/^\s*charset=(?:utf-8|us-ascii)\s*$/i.test(field))
    ) {
      throw new Error('Inline attachment declared media type has unsupported parameters.');
    }
    const embeddedCharset = fields.find((field) => /^charset=/i.test(field))?.toLowerCase();
    const declaredCharset = declaredFields[1]?.trim().toLowerCase();
    if (charset === undefined && declaredCharset !== undefined) charset = declaredCharset.slice(8);
    if (embeddedCharset && declaredCharset && embeddedCharset !== declaredCharset) {
      throw new Error('Inline attachment declared charset does not match its data URI.');
    }
  }
  const body = reference.uri.slice(comma + 1);
  const bytes = base64 ? decodeBase64(body) : decodePercent(body);
  if (charset === 'us-ascii' && bytes.some((byte) => byte > 127)) {
    throw new Error('Inline attachment bytes do not match declared us-ascii charset.');
  }
  if (reference.byteLength !== undefined && reference.byteLength !== bytes.length) {
    throw new Error('Inline attachment byteLength does not match its decoded bytes.');
  }
  const digest = await bytesDigest(bytes);
  if (reference.digest !== undefined && reference.digest !== digest) {
    throw new Error('Inline attachment digest does not match its decoded bytes (expected sha256:<hex>).');
  }
  return {
    mediaType: embeddedType,
    data: encodeBase64(bytes),
    byteLength: bytes.length,
    digest,
    ...(reference.name === undefined ? {} : { name: reference.name }),
  };
}

const textTypes = new Set([
  'text/plain',
  'text/markdown',
  'text/md',
  'text/html',
  'text/css',
  'text/xml',
  'text/csv',
  'text/rtf',
  'text/javascript',
  'application/json',
  'application/xml',
]);

export type VisibleProjector = (
  parts: readonly ModelVisiblePart[],
  model: string,
  placement: AttachmentPlacement,
) => Promise<readonly JsonObject[]>;

/** The `UnrepresentableInput` every provider's attachment projection throws. */
export function unsupportedAttachment(
  provider: string,
  model: string,
  type: string,
  placement: string,
): UnrepresentableInput {
  return new UnrepresentableInput(
    `${provider} model ${model} does not support inline ${type} attachments in ${placement} content.`,
  );
}

/** Text attachments have an explicit lossless UTF-8 mapping; binary types never fall back to text. */
export async function projectVisibleContent(
  parts: readonly ModelVisiblePart[],
  model: string,
  placement: AttachmentPlacement,
  project: (attachment: InlineAttachment, model: string, placement: AttachmentPlacement, text?: string) => JsonObject,
): Promise<readonly JsonObject[]> {
  const out: JsonObject[] = [];
  for (const part of parts) {
    if (part.kind === 'text') {
      out.push({ type: 'text', text: part.text });
      continue;
    }
    const attachment = await inlineAttachment(part.reference);
    if (attachment === undefined) {
      out.push({ type: 'text', text: `[ref ${part.reference.uri}]` });
      continue;
    }
    let text: string | undefined;
    if (textTypes.has(attachment.mediaType)) {
      try {
        text = strictUtf8.decode(decodeBase64(attachment.data));
      } catch {
        throw new Error('Inline text attachment must contain valid UTF-8 bytes.');
      }
    }
    out.push(project(attachment, model, placement, text));
  }
  return out;
}

/** Offline modality heuristics with a fixed margin and a transport-byte reserve. These rates
 * avoid charging encoded prompt characters; bytes cannot bound page/pixel/duration cost. */
export function estimateAttachmentTokens(part: JsonObject): number {
  if (typeof part.text === 'string') return 0; // Text is measured by the shared estimator.
  const bytes = typeof part.byteLength === 'number' ? part.byteLength : 0;
  const mime = typeof part.mediaType === 'string' ? part.mediaType : '';
  if (mime.startsWith('image/')) return 4096 + Math.ceil(bytes / 1024);
  if (mime === 'application/pdf') return 8192 + Math.ceil(bytes / 4);
  return 16384 + Math.ceil(bytes / 32);
}
