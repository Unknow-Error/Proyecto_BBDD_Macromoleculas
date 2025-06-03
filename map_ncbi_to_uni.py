def map_ncbi_to_uni(ncbi_id):
    run_url = "https://rest.uniprot.org/idmapping/run"
    params = {
        "from": "RefSeq_Protein",
        "to": "UniProtKB",
        "ids": ncbi_id
    }

    try:
        response = requests.post(run_url, data=params, timeout=30)
        response.raise_for_status()
        job_id = response.json()["jobId"]
    except Exception as e:
        raise RuntimeError(f"Error al enviar la solicitud de mapeo: {e}")

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
                    raise RuntimeError("El mapeo falló.")
            elif "results" in status_data:
                break
        except Exception as e:
            raise RuntimeError(f"Error al verificar el estado del trabajo: {e}")

    result_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"

    try:
        result_response = requests.get(result_url, timeout=30)
        result_response.raise_for_status()
        results_data = result_response.json()
    except Exception as e:
        raise RuntimeError(f"Error al obtener los resultados: {e}")

    uniprot_ids = [item["to"] for item in results_data.get("results", [])]

    return uniprot_ids
