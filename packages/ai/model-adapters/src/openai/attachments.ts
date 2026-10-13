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

const images = new Set(['image/jpeg', 'image/png', 'image/gif', 'image/webp']);
// Document formats accepted by Responses, mapped only through the documented native input_file
// path. Filename defaults keep the extension available to native document processing.
const openaiDocuments: Readonly<Record<string, string>> = {
  'application/pdf': 'pdf',
  'application/msword': 'doc',
  'application/vnd.openxmlformats-officedocument.wordprocessingml.document': 'docx',
  'application/rtf': 'rtf',
  'application/vnd.oasis.opendocument.text': 'odt',
  'application/vnd.ms-powerpoint': 'ppt',
  'application/vnd.openxmlformats-officedocument.presentationml.presentation': 'pptx',
  'application/vnd.ms-excel': 'xls',
  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet': 'xlsx',
};

export function projectVisibleContent(
  parts: readonly ModelVisiblePart[],
  model: string,
  placement: AttachmentPlacement,
): Promise<readonly JsonObject[]> {
  return projectCommonContent(parts, model, placement, (attachment, model, placement, text): JsonObject => {
    if (text !== undefined) return { type: 'text', text, attachment: { ...attachment } };
    const extension = openaiDocuments[attachment.mediaType];
    if (!((placement === 'user' || placement === 'tool-result') && (images.has(attachment.mediaType) || extension))) {
      throw unsupportedAttachment('openai', model, attachment.mediaType, placement);
    }
    return {
      type: 'file',
      ...attachment,
      ...(attachment.name === undefined && extension !== undefined ? { name: `attachment.${extension}` } : {}),
    };
  });
}
