/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { debounceTime } from 'rxjs/operators';
import { PluginStateAnimation } from '@molstar/plugin/state/animation/model';
import { PluginComponent } from '@molstar/plugin/state/component';
import { PluginContext } from '@molstar/plugin/context';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { encodeVideoAnimation } from '@molstar/video-export-extension/encoder';
import { MediabunnyBackend } from './backends/mediabunny.js';
import type { VideoRendererBackend } from './backend.js';

export interface VideoAnimationInfo {
  width: number;
  height: number;
}

function fpsParam() {
  return PD.Numeric(
    30,
    { min: 1, max: 120, step: 1 },
    {
      label: 'Frame rate (FPS)',
      help: () => ({
        description:
          'Number of frames rendered per second of animation. Higher values make motion smoother but require more rendering and encoding work. The animation duration stays the same. The default is 30 FPS.',
      }),
    },
  );
}

export class VideoControls extends PluginComponent {
  private paramsVersion = 0;
  private schemaKey = '';
  private disposed = false;
  private currentNames = new Set<string>();
  private animations: PluginStateAnimation[] = [];

  readonly behaviors = {
    animations: this.ev.behavior<PD.Params>({}),
    current: this.ev.behavior<{ anim: PluginStateAnimation; params: PD.Params; values: any } | undefined>(void 0),
    canApply: this.ev.behavior<PluginStateAnimation.CanApply>({ canApply: false }),
    info: this.ev.behavior<VideoAnimationInfo>({ width: 0, height: 0 }),
    params: this.ev.behavior<Record<string, any>>({ fps: 30 }),
    paramsSchema: this.ev.behavior<PD.Params>({ fps: fpsParam() }),
    encoder: this.ev.behavior<{ canApply: boolean; reason?: string }>({
      canApply: false,
      reason: 'Checking video codecs...',
    }),
  };

  setCurrent(name?: string) {
    const anim = this.animations.find((a) => a.name === name);
    if (!anim) {
      this.behaviors.current.next(anim);
      return;
    }

    const params = anim.params(this.plugin) as PD.Params;
    const values = PD.getDefaultValues(params);

    this.behaviors.current.next({ anim, params, values });
    this.behaviors.canApply.next(anim.canApply?.(this.plugin) ?? { canApply: true });
  }

  setCurrentParams(values: any) {
    this.behaviors.current.next({ ...this.behaviors.current.value!, values });
  }

  get current() {
    return this.behaviors.current.value;
  }

  render() {
    if (!this.behaviors.encoder.value.canApply) throw new Error(this.behaviors.encoder.value.reason);
    const task = Task.create('Export Animation', async (ctx) => {
      try {
        const resolution = this.plugin.helpers.viewportScreenshot?.getSizeAndViewport()!;
        const anim = this.current!;
        const movie = await encodeVideoAnimation(this.plugin, ctx, {
          animation: {
            definition: anim.anim,
            params: anim.values,
          },
          ...resolution,
          backend: this.backend,
          backendParams: this.behaviors.params.value,
          fps: this.behaviors.params.value.fps,
          pass: this.plugin.helpers.viewportScreenshot?.imagePass!,
        });

        const filename = anim.anim.display.name
          .toLowerCase()
          .replace(/\s/g, '-')
          .replace(/[^a-z0-9_\-]/g, '');
        return {
          movie: movie.data,
          mimeType: movie.mimeType,
          filename: `${this.plugin.helpers.viewportScreenshot?.getFilename('')}_${filename}.${movie.extension}`,
        };
      } catch (e) {
        this.plugin.log.error('Error during animation export');
        throw e;
      }
    });

    return this.plugin.runTask(task, { useOverlay: true });
  }

  private get manager() {
    return this.plugin.managers.animation;
  }

  private syncInfo() {
    const helper = this.plugin.helpers.viewportScreenshot;
    const size = helper?.getSizeAndViewport();
    if (!size) return;

    const alignment = this.backend.dimensionAlignment;
    this.behaviors.info.next({
      width: Math.floor(size.viewport.width / alignment) * alignment,
      height: Math.floor(size.viewport.height / alignment) * alignment,
    });
    void this.refreshParams();
  }

  private sync() {
    const animations = this.manager.animations.filter((a) => a.isExportable);

    const hasAll = animations.every((a) => this.currentNames.has(a.name));
    if (hasAll && this.currentNames.size === animations.length) {
      return;
    }

    const params = {
      current: PD.Select(
        animations[0]?.name,
        animations.map((a) => [a.name, a.display.name] as [string, string]),
        { label: 'Animation' },
      ),
    };

    const current = this.behaviors.current.value;
    const hasCurrent = !!animations.find((a) => a.name === current?.anim.name);

    this.currentNames = new Set(animations.map((a) => a.name));
    this.animations = animations;
    if (!hasCurrent) {
      this.setCurrent(animations[0]?.name);
    }
    this.behaviors.animations.next(params);
  }

  private init() {
    if (!this.plugin.canvas3d) return;

    this.subscribe(this.plugin.managers.animation.events.updated.pipe(debounceTime(16)), () => {
      this.sync();
    });

    this.subscribe(this.plugin.canvas3d.resized, () => this.syncInfo());
    this.subscribe(this.plugin.helpers.viewportScreenshot?.events.previewed!, () => this.syncInfo());

    this.subscribe(this.plugin.behaviors.state.isBusy, (b) => this.updateCanApply(b));
    this.subscribe(this.plugin.managers.snapshot.events.changed, (b) => this.updateCanApply(b));

    this.sync();
    this.syncInfo();
  }

  private updateCanApply(b?: any) {
    const anim = this.current;
    if (!b && anim) {
      this.behaviors.canApply.next(anim.anim.canApply?.(this.plugin) ?? { canApply: true });
    }
  }

  get videoParams() {
    const params = this.behaviors.paramsSchema.value;
    return PD.Group(params, {
      label: 'Video settings',
      isExpanded: true,
      presets: this.schemaKey ? this.backend.getPresets?.(params, this.behaviors.params.value) : undefined,
    });
  }

  setParams(values: Record<string, any>) {
    this.behaviors.params.next(values);
    void this.refreshParams();
  }

  private async refreshParams() {
    const version = ++this.paramsVersion;
    this.behaviors.encoder.next({ canApply: false, reason: 'Checking video codecs...' });
    const config = { ...this.behaviors.info.value, fps: this.behaviors.params.value.fps };
    try {
      if (config.width <= 0 || config.height <= 0) throw new Error('Video dimensions must be positive.');
      if (!Number.isFinite(config.fps) || config.fps <= 0) throw new Error('Frame rate must be positive.');
      const input = this.behaviors.params.value;
      const key = JSON.stringify([config.width, config.height, config.fps, input.container]);
      const params: PD.Params =
        key === this.schemaKey
          ? this.behaviors.paramsSchema.value
          : {
              fps: fpsParam(),
              ...(await this.backend.getParams(config, input)),
            };
      if (this.disposed || version !== this.paramsVersion) return;
      const values = PD.normalizeParams(params, input, 'all');
      this.schemaKey = JSON.stringify([config.width, config.height, config.fps, values.container]);
      this.behaviors.paramsSchema.next(params);
      this.behaviors.params.next(values);
      const reason = await this.backend.validate(config, values);
      if (this.disposed || version !== this.paramsVersion) return;
      this.behaviors.encoder.next({ canApply: !reason, reason });
    } catch (e) {
      if (this.disposed || version !== this.paramsVersion) return;
      this.behaviors.encoder.next({ canApply: false, reason: e instanceof Error ? e.message : String(e) });
    }
  }

  dispose() {
    this.disposed = true;
    ++this.paramsVersion;
    super.dispose();
  }

  constructor(
    private plugin: PluginContext,
    private backend: VideoRendererBackend = MediabunnyBackend,
  ) {
    super();
    this.init();
  }
}
