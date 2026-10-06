/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MarkdownExtension } from '../manager/markdown-extensions.js';
import {
  DisposeAudioMarkdownExtension,
  PauseAudioMarkdownExtension,
  PlayAudioMarkdownExtension,
  StopAudioMarkdownExtension,
  ToggleAudioMarkdownExtension,
} from './audio.js';
import { CenterCameraMarkdownExtension, FocusRefsMarkdownExtension } from './camera.js';
import { HighlightRefsMarkdownExtension } from './highlight.js';
import { QueryMarkdownExtension } from './query.js';
import {
  ApplySnapshotMarkdownExtension,
  NextSnapshotMarkdownExtension,
  PlaySnapshotsMarkdownExtension,
  PlayTransitionMarkdownExtension,
  StopAnimationMarkdownExtension,
} from './snapshots.js';

export const BuiltInMarkdownExtension: MarkdownExtension[] = [
  CenterCameraMarkdownExtension,
  ApplySnapshotMarkdownExtension,
  NextSnapshotMarkdownExtension,
  FocusRefsMarkdownExtension,
  HighlightRefsMarkdownExtension,
  QueryMarkdownExtension,
  PlayAudioMarkdownExtension,
  ToggleAudioMarkdownExtension,
  PauseAudioMarkdownExtension,
  StopAudioMarkdownExtension,
  DisposeAudioMarkdownExtension,
  PlayTransitionMarkdownExtension,
  PlaySnapshotsMarkdownExtension,
  StopAnimationMarkdownExtension,
];
