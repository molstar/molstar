/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useSeqViewer } from "@molstar/seq-viewer/react";
import { useSyncExternalStore } from "react";
import { TextAction } from "../ui";
import f from "./seqstar-format.module.css";
import { IssueList } from "./seqstar-issues";
import { patchNotices } from "./seqstar-patch";
import b from "./seqstar-patch.module.css";

/** The issues of the story's patches (`storyPatch`), under the viewer, each patch's own and dismissible. Inside the viewer's provider. */
export function PatchBanner() {
  const { context } = useSeqViewer();
  const store = patchNotices(context);
  const list = useSyncExternalStore(store.subscribe, () => store.list);
  if (list.length === 0) return null;
  return (
    <div className={b.banner} role="alert">
      {list.map((n) => (
        <div key={n.id} className={f.stack}>
          <div className={f.row}>
            <span className={f.code}>"{n.label}" came back with issues, so it may not have done what it says:</span>
            <TextAction tone="quiet" onClick={() => store.dismiss(n.id)}>
              Dismiss
            </TextAction>
          </div>
          <IssueList issues={n.issues} />
        </div>
      ))}
    </div>
  );
}
