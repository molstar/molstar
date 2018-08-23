#if defined(dColorType_attribute)
    vColor.rgb = aColor;
#elif defined(dColorType_instance)
    vColor.rgb = readFromTexture(tColor, aInstance, uColorTexDim).rgb;
#elif defined(dColorType_group)
    vColor.rgb = readFromTexture(tColor, aGroup, uColorTexDim).rgb;
#elif defined(dColorType_groupInstance)
    vColor.rgb = readFromTexture(tColor, aGroup * float(uGroupCount) + aGroup, uColorTexDim).rgb;
#elif defined(dColorType_objectPicking)
    vColor = encodeIdRGBA(float(uObjectId));
#elif defined(dColorType_instancePicking)
    vColor = encodeIdRGBA(aInstance);
#elif defined(dColorType_groupPicking)
    vColor = encodeIdRGBA(aGroup);
#endif