import React from 'react';
import { AnimateCameraSpin } from '../../mol-plugin-state/animation/built-in/camera-spin';
import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { Button } from '../../mol-plugin-ui/controls/common';
import { Task } from '../../mol-task';
import { download } from '../../mol-util/download';
import { encodeMp4Animation, Mp4Encoder } from './encoder';

interface State {
    data?: { movie: Uint8Array };
}

export class Mp4EncoderTestUI extends CollapsableControls<{}, State> {
    protected defaultState(): State & CollapsableState  {
        return {
            header: 'Export MP4',
            isCollapsed: false,
            brand: { accent: 'cyan' }
        };
    }
    protected renderControls(): JSX.Element | null {
        return <>
            <Button onClick={() => this.gen()}>Generate</Button>
            {this.state.data && <Button onClick={() => this.save()}>Save</Button>}
        </>;
    }

    save() {
        download(new Blob([this.state.data!.movie]), 'test.mp4');
        // download(this.state.data!.image, 'test.png');
    }

    gen() {
        const task = Task.create('Export Animation', async ctx => {
            const movie = await encodeMp4Animation(this.plugin, ctx, {
                animation: {
                    definition: AnimateCameraSpin,
                    params: { durationInMs: 2000, speed: 1, direction: 'cw' }
                },
                width: 1280,
                height: 720,
                quantizationParameter: 18,
                pass: this.plugin.helpers.viewportScreenshot?.imagePass!,
            });

            this.setState({ data: { movie } });
        });

        this.plugin.runTask(task);
    }

    async generate() {
        const encoder = new Mp4Encoder(this.plugin);
        const data = await encoder.generate();
        this.setState({ data });
    }
}