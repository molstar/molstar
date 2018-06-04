#if defined(dColorType_attribute)
    vColor.rgb = aColor;
#elif defined(dColorType_instance)
    vColor.rgb = readFromTexture(tColor, aInstanceId, uColorTexSize).rgb;
#elif defined(dColorType_element)
    vColor.rgb = readFromTexture(tColor, aElementId, uColorTexSize).rgb;
#elif defined(dColorType_elementInstance)
    vColor.rgb = readFromTexture(tColor, aInstanceId * float(uElementCount) + aElementId, uColorTexSize).rgb;
#elif defined(dColorType_objectPicking)
    vColor = encodeIdRGBA(float(uObjectId));
#elif defined(dColorType_instancePicking)
    vColor = encodeIdRGBA(aInstanceId);
#elif defined(dColorType_elementPicking)
    vColor = encodeIdRGBA(aElementId);
#endif