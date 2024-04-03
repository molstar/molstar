import React, { useEffect, useRef } from 'react';
import JSONEditor, { JSONEditorOptions } from 'jsoneditor';
import 'jsoneditor/dist/jsoneditor.css';
import { Card } from 'primereact/card';
import { VolsegEntryData } from './entry-root';
import { Button } from '../../../mol-plugin-ui/controls/common';
import { AnnotationMetadata } from './volseg-api/data';

interface JSONEditorComponentProps {
    jsonData: any;
    entryData: VolsegEntryData
}

async function updateJSON(jsonData: AnnotationMetadata, entryData: VolsegEntryData) {
    console.log('JSON changed'); console.log(jsonData);
    await entryData.api.updateAnnotationsJson(entryData.source, entryData.entryId, jsonData);
    await entryData.updateMetadata();
}

export const JSONEditorComponent: React.FC<JSONEditorComponentProps> = ({ jsonData, entryData }) => {
    const containerRef = useRef<HTMLDivElement | null>(null);
    let jsonEditor: JSONEditor | null = null;
    // let jsonDataUpdated = jsonData;
    const jsonDataUpdated = useRef(jsonData);
    useEffect(() => {
        if (containerRef.current) {
            const options: JSONEditorOptions = {
                onChangeJSON: async (jsonData) => {
                    jsonDataUpdated.current = jsonData;
                }
            };
            jsonEditor = new JSONEditor(containerRef.current, options);
            jsonEditor.set(jsonData);
        }

        return () => {
            if (jsonEditor) {
                jsonEditor.destroy();
            }
        };
    }, [jsonData]);


    return (
        <Card>
            <div ref={containerRef} style={{ width: '100%', height: '400px' }} />
            <Button onClick={async () => await updateJSON(jsonDataUpdated.current, entryData)}>Update JSON</Button>
        </Card>
    );
};