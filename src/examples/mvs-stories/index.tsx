/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getMolComponentContext } from './context';
import './index.html';
import './elements/snapshot-markdown';
import './elements/viewer';
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

type Story = { kind: 'built-in', id: string } | { kind: 'url', url: string, format: 'mvsx' | 'mvsj' } | undefined;
const CurrentStory = new BehaviorSubject<Story>(undefined);

function SelectStoryUI({ subject }: { subject: BehaviorSubject<Story> }) {
    const current = useBehavior(subject);
    const selectedId = current?.kind === 'built-in' ? current.id : current?.kind === 'url' ? 'url' : '';

    return <select onChange={e => {
        const value = e.currentTarget.value;
        const s = Stories.find(s => s.id === value);
        if (!s) return;
        subject.next({ kind: 'built-in', id: s.id });
    }} value={selectedId}>
        {!current && <option value=''>Select a story...</option>}
        {Stories.map(s => <option key={s.name} value={s.id}>Story: {s.name}</option>)}
        {current?.kind === 'url' && <option disabled>------------------</option>}
        {current?.kind === 'url' && <option value='url'>{current.url}</option>}
    </select>;
}

function init() {
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

    createRoot(document.getElementById('select-story')!).render(<SelectStoryUI subject={CurrentStory} />);
}

(window as any).mc = MC;
(window as any).downloadStory = () => {
    if (CurrentStory.value?.kind !== 'built-in') return;
    const id = CurrentStory.value.id;
    const story = Stories.find(s => s.id === id);
    if (!story) return;
    const data = JSON.stringify(story.buildStory(), null, 2);
    download(new Blob([data], { type: 'application/json' }), 'story.mvsj');
};
(window as any).initStories = init;
(window as any).CurrentStory = CurrentStory;