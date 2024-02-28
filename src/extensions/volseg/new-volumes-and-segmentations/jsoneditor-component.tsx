import React, { useEffect, useRef } from "react";
import JSONEditor, { JSONEditorOptions } from "jsoneditor";
import "jsoneditor/dist/jsoneditor.css";
import { Card } from "primereact/card";
import { VolsegEntryData } from "./entry-root";
import { MetadataWrapper } from "./volseg-api/utils";
import { Button } from "../../../mol-plugin-ui/controls/common";

interface JSONEditorComponentProps {
    jsonData: any;
    entryData: VolsegEntryData
}

async function updateJSON(jsonData, entryData) {
    console.log('JSON changed'); console.log(jsonData);
    await entryData.api.updateAnnotationsJson(entryData.source, entryData.entryId, jsonData);
    // should fetch metadata again
    // const metadata = await entryData.api.getMetadata(entryData.source, entryData.entryId);
    // entryData.metadata.value! = new MetadataWrapper(metadata);
    debugger;
    await entryData.updateMetadata();
}

const JSONEditorComponent: React.FC<JSONEditorComponentProps> = ({ jsonData, entryData }) => {
    const containerRef = useRef<HTMLDivElement | null>(null);
    let jsonEditor: JSONEditor | null = null;
    // let jsonDataUpdated = jsonData;
    const jsonDataUpdated = useRef(jsonData);
    useEffect(() => {
        if (containerRef.current) {
            const options: JSONEditorOptions = {
                onChangeJSON: async (jsonData) => {
                    // console.log('JSON changed'); console.log(jsonData);
                    // await entryData.api.updateAnnotationsJson(entryData.source, entryData.entryId, jsonData);
                    // // should fetch metadata again
                    // // const metadata = await entryData.api.getMetadata(entryData.source, entryData.entryId);
                    // // entryData.metadata.value! = new MetadataWrapper(metadata);
                    
                    // await entryData.updateMetadata();
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


    // TODO: can try to use .get() method
    // const jsonContent = containerRef.current,options.get()
    return (
        <Card>
            <div ref={containerRef} style={{ width: '100%', height: '400px' }} />
            <Button onClick={async () => await updateJSON(jsonDataUpdated.current, entryData)}>Update JSON</Button>
        </Card>
    );
};

export default JSONEditorComponent;