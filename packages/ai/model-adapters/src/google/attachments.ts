/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { JsonObject, ModelVisiblePart } from '@molstar/ai-agent';
import {
  type AttachmentPlacement,
  projectVisibleContent as projectCommonContent,
  unsupportedAttachment,
} from '../common/attachments.js';

const googleMedia = new Set([
  'application/pdf',
  'image/png',
  'image/jpeg',
  'image/webp',
  'image/heic',
  'image/heif',
  'audio/wav',
  'audio/mp3',
  'audio/mpeg',
  'audio/aiff',
  'audio/aac',
  'audio/ogg',
  'audio/flac',
  'video/mp4',
  'video/mpeg',
  'video/quicktime',
  'video/avi',
  'video/x-flv',
  'video/mpg',
  'video/webm',
  'video/wmv',
  'video/3gpp',
]);

/** Gemini 3 and later (`gemini-3`, `gemini-3.8-flash`, ...). */
export const isGemini3 = (model: string): boolean => /^gemini-(?:[3-9]|[1-9]\d)(?:[.-]|$)/.test(model);

export function projectVisibleContent(
  parts: readonly ModelVisiblePart[],
  model: string,
  placement: AttachmentPlacement,
): Promise<readonly JsonObject[]> {
  return projectCommonContent(parts, model, placement, (attachment, model, placement, text): JsonObject => {
    if (text !== undefined) return { type: 'text', text, attachment: { ...attachment } };
    if (
      !(
        placement !== 'system' &&
        googleMedia.has(attachment.mediaType) &&
        (placement !== 'tool-result' || isGemini3(model))
      )
    ) {
      throw unsupportedAttachment('google', model, attachment.mediaType, placement);
    }
    return { type: 'file', ...attachment };
  });
}
