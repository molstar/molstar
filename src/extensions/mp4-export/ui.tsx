/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import React from 'react';
import { AnimateCameraSpin } from '../../mol-plugin-state/animation/built-in/camera-spin';
import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { Button } from '../../mol-plugin-ui/controls/common';
import { SubscriptionsOutlinedSvg } from '../../mol-plugin-ui/controls/icons';
import { Task } from '../../mol-task';
import { download } from '../../mol-util/download';
import { encodeMp4Animation, Mp4Encoder } from './encoder';

interface State {
    data?: { movie: Uint8Array };
}

export class Mp4EncoderUI extends CollapsableControls<{}, State> {
    protected defaultState(): State & CollapsableState  {
        return {
            header: 'Export Animation',
            isCollapsed: true,
            brand: { accent: 'cyan', svg: SubscriptionsOutlinedSvg }
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
    }

    gen() {
        const task = Task.create('Export Animation', async ctx => {
            const resolution = this.plugin.helpers.viewportScreenshot?.getSizeAndViewport()!;
            const movie = await encodeMp4Animation(this.plugin, ctx, {
                animation: {
                    definition: AnimateCameraSpin,
                    params: { durationInMs: 2000, speed: 1, direction: 'cw' }
                },
                ...resolution,
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