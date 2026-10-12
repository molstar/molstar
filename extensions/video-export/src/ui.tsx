/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { merge } from 'rxjs';
import { debounceTime } from 'rxjs/operators';
import { CollapsableControls, type CollapsableState } from '@molstar/plugin-ui/controls/collapsable';
import { Button } from '@molstar/plugin-ui/controls/common';
import { CameraOutlinedSvg, GetAppSvg, Icon, SubscriptionsOutlinedSvg } from '@molstar/plugin-ui/controls/icons';
import { ParameterControls } from '@molstar/plugin-ui/controls/parameters';
import { download } from '@molstar/core/util/download';
import { VideoControls } from '@molstar/video-export-extension/controls';

interface State {
  busy?: boolean;
  data?: { movie: Uint8Array<ArrayBuffer>; filename: string; mimeType: string };
}

export class VideoEncoderUI extends CollapsableControls<{}, State> {
  private _controls: VideoControls | undefined;

  get controls() {
    return this._controls || (this._controls = new VideoControls(this.plugin));
  }

  protected defaultState(): State & CollapsableState {
    return {
      header: 'Export Animation',
      isCollapsed: true,
      brand: { accent: 'cyan', svg: SubscriptionsOutlinedSvg },
    };
  }

  private downloadControls() {
    return (
      <>
        <div className="msp-control-offset msp-help-text">
          <div className="msp-help-description" style={{ textAlign: 'center' }}>
            Rendering successful!
          </div>
        </div>
        <Button icon={GetAppSvg} onClick={this.save} style={{ marginTop: 1 }}>
          Save Animation
        </Button>
        <Button onClick={() => this.setState({ data: void 0 })} style={{ marginTop: 6 }}>
          Clear
        </Button>
      </>
    );
  }

  protected renderControls(): JSX.Element | null {
    if (this.state.data) {
      return this.downloadControls();
    }

    const ctrl = this.controls;
    const current = ctrl.behaviors.current.value;
    const info = ctrl.behaviors.info.value;
    const animation = ctrl.behaviors.canApply.value;
    const encoder = ctrl.behaviors.encoder.value;
    const canApply = animation.canApply ? encoder : animation;
    return (
      <>
        <ParameterControls
          params={ctrl.behaviors.animations.value}
          values={{ current: current?.anim.name }}
          onChangeValues={(xs) => ctrl.setCurrent(xs.current)}
          isDisabled={this.state.busy}
        />
        {current && (
          <ParameterControls
            params={current.params}
            values={current.values}
            onChangeValues={(xs) => ctrl.setCurrentParams(xs)}
            isDisabled={this.state.busy}
          />
        )}
        <div className="msp-control-offset msp-help-text">
          <div className="msp-help-description" style={{ textAlign: 'center' }}>
            Resolution: {info.width}x{info.height}
            <br />
            Adjust in viewport using <Icon svg={CameraOutlinedSvg} inline />
          </div>
        </div>
        <ParameterControls
          params={{ video: ctrl.videoParams }}
          values={{ video: ctrl.behaviors.params.value }}
          onChangeValues={(xs) => ctrl.setParams(xs.video)}
          isDisabled={this.state.busy}
        />
        <Button
          onClick={this.generate}
          style={{ marginTop: 1 }}
          disabled={this.state.busy || !canApply.canApply}
          commit={canApply.canApply ? 'on' : 'off'}
        >
          {canApply.canApply ? 'Render' : (canApply.reason ?? 'Invalid params/state')}
        </Button>
      </>
    );
  }

  componentDidMount() {
    const merged = merge(
      this.controls.behaviors.animations,
      this.controls.behaviors.current,
      this.controls.behaviors.canApply,
      this.controls.behaviors.info,
      this.controls.behaviors.params,
      this.controls.behaviors.paramsSchema,
      this.controls.behaviors.encoder,
    );

    this.subscribe(merged.pipe(debounceTime(10)), () => {
      if (!this.state.isCollapsed) this.forceUpdate();
    });
  }

  componentWillUnmount() {
    super.componentWillUnmount();
    this._controls?.dispose();
    this._controls = void 0;
  }

  save = () => {
    download(new Blob([this.state.data!.movie], { type: this.state.data!.mimeType }), this.state.data!.filename);
  };

  generate = async () => {
    try {
      this.setState({ busy: true });
      const data = await this.controls.render();
      this.setState({ busy: false, data });
    } catch (e) {
      console.error(e);
      this.setState({ busy: false });
    }
  };
}
