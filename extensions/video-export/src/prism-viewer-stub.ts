/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LazyViewerExtension, ViewerSkill } from '@molstar/prism-viewer/spec';

// The video export as a lazy Prism extension (viewer-spec.md §3, parity.md §6: Export › Animation; lazy): this module
// is the stub, plain data and types only, so listing the extension pulls in neither the encoder nor the views. Its
// operations and commands exist from the start; the first use loads `prism-viewer.ts`, checked against this.

/** The operations, with the flags their definitions have (the host compares them when the extension loads). */
export const VideoOperations = {
  export: {
    id: 'video.export',
    title: 'Export a video',
    history: false,
    undoable: false,
    modules: ['canvas'],
  },
  cancel: { id: 'video.cancel', title: 'Cancel the video export', history: false },
} as const;

export const VideoSection = 'video-export.section.animation';

export const VideoCommands = {
  export: 'video.export',
  cancel: 'video.cancel',
  open: 'video.open',
} as const;

export const VideoSkill: ViewerSkill = {
  id: 'video',
  title: 'Video export',
  summary: 'record an animation (spin, rock, snapshots) to an MP4 or WebM file',
  guidance: [
    'Skill video: export_video renders an animation frame by frame into a video file (MP4 or WebM).',
    '- animation: a registered name (animation.get lists them): built-in.animate-camera-spin,',
    '  built-in.animate-camera-rock, built-in.animate-state-snapshots, …; params are its own.',
    "- duration (seconds) is the video's length (spin and rock play over it); leave it out for the animation's own.",
    '- It needs a canvas, and takes as long as the frames take to render; cancel_video_export stops it.',
  ].join('\n'),
  tools: [
    { id: 'export_video', operation: 'video.export', description: 'Records an animation to a video file.' },
    { id: 'cancel_video_export', operation: 'video.cancel', description: 'Stops a video export in progress.' },
  ],
};

export const LazyVideoExport: LazyViewerExtension = {
  id: 'video-export',
  title: 'Video export',
  stub: {
    parts: [
      { id: VideoSection, kind: 'section', slot: 'export.picture', title: 'Animation', icon: 'clapperboard', order: 2 },
    ],
    operations: [VideoOperations.export, VideoOperations.cancel],
    commands: [
      { id: VideoCommands.export, title: 'Export a video…', category: 'Export', icon: 'clapperboard' },
      { id: VideoCommands.cancel, title: 'Cancel the video export', category: 'Export', icon: 'x' },
      { id: VideoCommands.open, title: 'Show the Animation export', category: 'Export', icon: 'clapperboard' },
    ],
    skills: [VideoSkill],
  },
  load: async () => (await import('./prism-viewer.js')).VideoExportExtension,
};
