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

const anthropicImages = new Set(['image/jpeg', 'image/png', 'image/gif', 'image/webp']);

export function projectVisibleContent(
  parts: readonly ModelVisiblePart[],
  model: string,
  placement: AttachmentPlacement,
): Promise<readonly JsonObject[]> {
  return projectCommonContent(parts, model, placement, (attachment, model, placement, text): JsonObject => {
    if (text !== undefined) {
      if (placement === 'user')
        return { type: 'file', ...attachment, mediaType: 'text/plain', sourceMediaType: attachment.mediaType, text };
      return { type: 'text', text, attachment: { ...attachment } };
    }
    if (
      !(
        (placement === 'user' || placement === 'tool-result') &&
        (anthropicImages.has(attachment.mediaType) || attachment.mediaType === 'application/pdf')
      )
    ) {
      throw unsupportedAttachment('anthropic', model, attachment.mediaType, placement);
    }
    return { type: 'file', ...attachment };
  });
}
