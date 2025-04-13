/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getMolComponentContext } from './context';
import './index.html';
import './elements/snapshot-markdown';
import './elements/viewer';
import { buildStory } from './stories/kinase';
import '../../mol-plugin-ui/skin/light.scss';
import './styles.scss';
import { download } from '../../mol-util/download';
import { BehaviorSubject } from 'rxjs';
import { Stories } from './stories';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';
import { createRoot } from 'react-dom/client';

export class MolComponents {
    getContext(name?: string) {
        return getMolComponentContext({ name });
    }
}

const MC = new MolComponents();

const CurrentStory = new BehaviorSubject<{ kind: 'built-in', id: string } | { kind: 'url', url: string, format: 'mvsx' | 'mvsj' } | undefined>(undefined);
CurrentStory.subscribe(story => {
    if (!story) {
        history.replaceState({}, '', '');
    } else if (story.kind === 'url') {
        history.replaceState({}, '', story ? `?story-url=${encodeURIComponent(story.url)}&data-format=${story.format}` : '');
        MC.getContext().dispatch({
            kind: 'load-mvs',
            format: story.format,
            url: story.url,
        });
    } else if (story.kind === 'built-in') {
        history.replaceState({}, '', story ? `?story=${story.id}` : '');
        const s = Stories.find(s => s.id === story.id);
        if (s) {
            MC.getContext().dispatch({
                kind: 'load-mvs',
                data: s.buildStory(),
            });
        } else {
            console.warn('Story not found:', story.id);
            CurrentStory.next({ kind: 'built-in', id: Stories[0].id });
        }
    }
});

function SelectStoryUI() {
    const current = useBehavior(CurrentStory);

    return <select onChange={e => {
        const value = e.currentTarget.value;
        const s = Stories.find(s => s.id === value);
        if (!s) return;
        CurrentStory.next({ kind: 'built-in', id: s.id });
    }}>
        {!current && <option value=''>Select a story...</option>}
        {Stories.map(s => <option key={s.name} value={s.id} selected={current?.kind === 'built-in' && current.id === s.id}>Story: {s.name}</option>)}
        {current?.kind === 'url' && <option disabled>------------------</option>}
        {current?.kind === 'url' && <option value='url' selected>{current.url}</option>}
    </select>;
}

(window as any).mc = MC;
(window as any).buildStory = buildStory;
(window as any).downloadStory = () => {
    if (CurrentStory.value?.kind !== 'built-in') return;
    const name = CurrentStory.value.id;
    const story = Stories.find(s => s.name === name);
    const data = JSON.stringify(story, null, 2);
    download(new Blob([data], { type: 'application/json' }), 'story.mvsj');
};
(window as any).init = () => {
    const urlParams = new URLSearchParams(window.location.search);
    const storyUrl = urlParams.get('story-url');
    const dataFormat = urlParams.get('data-format') as 'mvsx' | 'mvsj' | null;
    const storyId = urlParams.get('story');

    if (storyUrl) {
        CurrentStory.next({ kind: 'url', url: storyUrl, format: dataFormat ?? 'mvsj' });
    } else if (storyId) {
        CurrentStory.next({ kind: 'built-in', id: storyId });
    } else {
        CurrentStory.next({ kind: 'built-in', id: Stories[0].id });
    }

    createRoot(document.getElementById('select-story')!).render(<SelectStoryUI />);
};