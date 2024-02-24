import React, { useEffect, useRef } from "react";
import JSONEditor, { JSONEditorOptions } from "jsoneditor";
import "jsoneditor/dist/jsoneditor.css";
import { Card } from "primereact/card";
import { VolsegEntryData } from "./entry-root";
import { MetadataWrapper } from "./volseg-api/utils";

interface JSONEditorComponentProps {
    jsonData: any;
    entryData: VolsegEntryData
}

const JSONEditorComponent: React.FC<JSONEditorComponentProps> = ({ jsonData, entryData }) => {
    const containerRef = useRef<HTMLDivElement | null>(null);
    let jsonEditor: JSONEditor | null = null;

    useEffect(() => {
        if (containerRef.current) {
            const options: JSONEditorOptions = {
                onChangeJSON: async (jsonData) => {
                    console.log('JSON changed'); console.log(jsonData);
                    await entryData.api.updateAnnotationsJson(entryData.source, entryData.entryId, jsonData);
                    // should fetch metadata again
                    // const metadata = await entryData.api.getMetadata(entryData.source, entryData.entryId);
                    // entryData.metadata = new MetadataWrapper(metadata);
                    
                    await entryData.updateMetadata();
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
        </Card>
    );
};

export default JSONEditorComponent;