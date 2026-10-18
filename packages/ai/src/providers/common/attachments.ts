/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ModelVisiblePart, Reference } from '../../messages/index.js';
import { hashing, UnrepresentableInput } from '../../messages/index.js';

// Inline attachments (`data:` URIs) as the AI SDK's file parts. Which media a provider takes where is its rule
// (`AttachmentRule`); what it cannot take is `UnrepresentableInput` at projection, so the harness can say so before any
// request. Corrupt bytes (length, digest, encoding) are plain errors: the input is wrong, not unrepresentable.

export type AttachmentPlacement = 'system' | 'user' | 'assistant' | 'tool-result';

/** Where a provider takes binary media, and in what form. */
export interface AttachmentRule {
  /** The name in the error ("openai-chat"). */
  readonly name: string;
  accepts(mediaType: string, placement: AttachmentPlacement, model: string): boolean;
  /** Text attachments as a text document (Anthropic's, in user content); else as text. */
  textAsFile?(placement: AttachmentPlacement): boolean;
  /** A file name for a type the provider needs one for (OpenAI's documents). */
  filename?(mediaType: string): string | undefined;
}

/** The parts a provider takes: text, or a file (base64 bytes or a text document). JSON, as the payload retains it. */
export type VisiblePart =
  | { readonly type: 'text'; readonly text: string }
  | {
      readonly type: 'file';
      readonly data:
        | { readonly type: 'data'; readonly data: string }
        | { readonly type: 'text'; readonly text: string };
      readonly mediaType: string;
      readonly filename?: string;
    };

export interface Visible {
  readonly parts: readonly VisiblePart[];
  /** Whether an attachment is among them: the parts then stay parts, else they are one text. */
  readonly rich: boolean;
  readonly text: string;
}

/** An inline attachment's bytes, decoded and checked. */
export interface InlineAttachment {
  readonly mediaType: string;
  readonly data: string;
  readonly byteLength: number;
  readonly name?: string;
}

const strictUtf8 = new TextDecoder('utf-8', { fatal: true, ignoreBOM: true });

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

/** `murmur3-128:<32 hex>` of inline bytes, the form a reference's `digest` takes (the messages' own hash). */
export function bytesDigest(bytes: Uint8Array): string {
  return `murmur3-128:${hashing.bytes(bytes)}`;
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

/** A `data:` URI's bytes, checked against what the reference says of them; undefined for any other URI. */
export function inlineAttachment(reference: Reference): InlineAttachment | undefined {
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
  if (reference.digest !== undefined && reference.digest !== bytesDigest(bytes)) {
    throw new Error('Inline attachment digest does not match its decoded bytes (expected murmur3-128:<hex>).');
  }
  return {
    mediaType: embeddedType,
    data: encodeBase64(bytes),
    byteLength: bytes.length,
    ...(reference.name === undefined ? {} : { name: reference.name }),
  };
}

/** Text and inline attachments as the provider takes them where they stand; other references as `[ref <uri>]`. */
export function visibleParts(
  parts: readonly ModelVisiblePart[],
  placement: AttachmentPlacement,
  model: string,
  rule: AttachmentRule,
): Visible {
  const out: VisiblePart[] = [];
  let rich = false;
  for (const part of parts) {
    if (part.kind === 'text') {
      out.push({ type: 'text', text: part.text });
      continue;
    }
    const attachment = inlineAttachment(part.reference);
    if (attachment === undefined) {
      out.push({ type: 'text', text: `[ref ${part.reference.uri}]` });
      continue;
    }
    rich = true;
    const filename = attachment.name ?? rule.filename?.(attachment.mediaType);
    if (textTypes.has(attachment.mediaType)) {
      let text: string;
      try {
        text = strictUtf8.decode(decodeBase64(attachment.data));
      } catch {
        throw new Error('Inline text attachment must contain valid UTF-8 bytes.');
      }
      out.push(
        rule.textAsFile?.(placement)
          ? {
              type: 'file',
              data: { type: 'text', text },
              mediaType: 'text/plain',
              ...(attachment.name === undefined ? {} : { filename: attachment.name }),
            }
          : { type: 'text', text },
      );
      continue;
    }
    if (!rule.accepts(attachment.mediaType, placement, model)) {
      throw new UnrepresentableInput(
        `${rule.name} model ${model} does not support inline ${attachment.mediaType} attachments in ${placement} content.`,
      );
    }
    out.push({
      type: 'file',
      data: { type: 'data', data: attachment.data },
      mediaType: attachment.mediaType,
      ...(filename === undefined ? {} : { filename }),
    });
  }
  return { parts: out, rich, text: out.map((part) => (part.type === 'text' ? part.text : '')).join('') };
}

/** Offline modality heuristics with a fixed margin and a transport-byte reserve. These rates
 * avoid charging encoded prompt characters; bytes cannot bound page/pixel/duration cost. */
export function estimateAttachmentTokens(mediaType: string, bytes: number): number {
  if (mediaType.startsWith('image/')) return 4096 + Math.ceil(bytes / 1024);
  if (mediaType === 'application/pdf') return 8192 + Math.ceil(bytes / 4);
  return 16384 + Math.ceil(bytes / 32);
}

/** The byte length of canonical base64. */
export function base64Bytes(data: string): number {
  return Math.floor((data.length * 3) / 4) - (data.endsWith('==') ? 2 : data.endsWith('=') ? 1 : 0);
}
