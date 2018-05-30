#if defined(dColorType_attribute)
    vColor.xyz = aColor;
#elif defined(dColorType_instance)
    vColor.xyz = read_vec3(tColor, aInstanceId, uColorTexSize);
#elif defined(dColorType_element)
    vColor.xyz = read_vec3(tColor, aElementId, uColorTexSize);
#elif defined(dColorType_elementInstance)
    vColor.xyz = read_vec3(tColor, aInstanceId * float(uElementCount) + aElementId, uColorTexSize);
#elif defined(dColorType_objectPicking)
    vColor = encodeIdRGBA(float(uObjectId));
#elif defined(dColorType_instancePicking)
    vColor = encodeIdRGBA(aInstanceId);
#elif defined(dColorType_elementPicking)
    vColor = encodeIdRGBA(aElementId);
#endif