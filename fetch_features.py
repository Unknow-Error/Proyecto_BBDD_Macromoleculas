def fetch_features(identifier):
    url = f"https://rest.uniprot.org/uniprotkb/{identifier}.json"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    protein_data = response.json()
    features = protein_data.get("features", [])

    feature_list = []
    for feature in features:
        loc = feature.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")

        if start and end:
            feature_info = {
                "type": feature.get("type", "Unknown"),
                "description": feature.get("description", "No description"),
                "start": int(start),
                "end": int(end)
            }
            feature_list.append(feature_info)
    return feature_list
