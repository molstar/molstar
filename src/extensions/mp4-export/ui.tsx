import React from 'react';
import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { Button } from '../../mol-plugin-ui/controls/common';
import { download } from '../../mol-util/download';
import { Mp4Encoder } from './encoder';

interface State {
    data?: { movie: Uint8Array, image: Blob };
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
            <Button onClick={() => this.generate()}>Generate</Button>
            {this.state.data && <Button onClick={() => this.save()}>Save</Button>}
        </>;
    }

    save() {
        download(new Blob([this.state.data!.movie]), 'test.mp4');
        // download(this.state.data!.image, 'test.png');
    }

    async generate() {
        const encoder = new Mp4Encoder(this.plugin);
        const data = await encoder.generate();
        this.setState({ data });
    }
}