def fetch_structures(identifier):
    base_url = "https://rest.uniprot.org/idmapping/run"
    params = {
        "from": "UniProtKB_AC-ID",
        "to": "PDB",
        "ids": identifier
    }

    response = requests.post(base_url, data=params, timeout=30)
    response.raise_for_status()
    job_id = response.json()["jobId"]

    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        try:
            time.sleep(3)
            status_response = requests.get(status_url, timeout=30)
            status_response.raise_for_status()
            status_data = status_response.json()

            if "jobStatus" in status_data:
                job_status = status_data["jobStatus"]
                if job_status == "FINISHED":
                    break
                elif job_status == "FAILED":
                    raise RuntimeError("El mapeo UniProt ↔ PDB falló..")
            elif "results" in status_data:
                break
        except Exception as e:
            raise RuntimeError(f"Error al verificar el estado del trabajo: {e}")

    results_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    results_response = requests.get(results_url, timeout=30)
    results_response.raise_for_status()
    results_data = results_response.json()

    pdb_ids = []
    for result in results_data.get("results", []):
        pdb_id = result.get("to")
        if pdb_id:
            pdb_ids.append(pdb_id)

    return pdb_ids
