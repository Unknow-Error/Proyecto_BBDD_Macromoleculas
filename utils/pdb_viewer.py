import webbrowser
import py3Dmol
import random
import re
from matplotlib.colors import CSS4_COLORS
import requests
import pandas as pd
from io import StringIO
import time
import os
import shutil

URL_Uniprot = "https://rest.uniprot.org/uniprotkb/search"

class API_Uniprot_rest:
  """
    Cliente para la REST API moderna de UniProt u otras como UniRef.
    Devuelve resultados en formato TSV como pandas.DataFrame.
  """
  def __init__(self, parametros, base_url, timeout=30, cursor=None):
    self.parametros = parametros
    self.base_url = base_url
    self.timeout = timeout

  def _conseguir_pagina(self, url, parametros):
    """
      Realiza la petición HTTP a UniProt  y retorna el DataFrame de la página y el link next (o None).
    """
    #Peticion:
    response = requests.get(url, params=parametros, timeout=self.timeout)
    response.raise_for_status()
    #Leer TSV:
    df = pd.read_csv(StringIO(response.content.decode("utf-8")), sep="\t", header=0)
    # Extraer enlace "siguiente" de encabezado Link
    siguiente_link = None
    link_header = response.headers.get("Link", "")
    # Updated regex pattern to extract the URL with the scheme
    pattern = re.compile(r'<([^>]+)>;\s*rel="next"')
    match = pattern.search(link_header)
    if match:
        siguiente_link = match.group(1)

    return df, siguiente_link

  def data_de_paginacion_tsv(self, delay):
    '''
      Similar a response_data pero paginable, utilizando la funcion anterior.
      Retorna un DataFrame con los datos obtenidos en batch de la paginacion.
      Requiere que en parametros se solicite los datos en .tsv
    '''
    # Preparar primera petición
    url_actual = self.base_url
    parametros = self.parametros.copy()

    paginas = []  # Para guardar de cada pagina

    while True:
      time.sleep(delay)
      df, siguiente_link = self._conseguir_pagina(url_actual, parametros)
      if df.empty:
        break
      paginas.append(df) # Guardar la página
      if siguiente_link is None:
        break
      # Preparar siguiente iteración: usar URL completa en next_link
      url_actual  = siguiente_link
      # Vaciar params para que requests use URL tal cual
      parametros = None

    # Concatenamos todo en un solo DataFrame (si hubo al menos una página)
    return pd.concat(paginas, ignore_index=True) if paginas else pd.DataFrame()

class PDB_Viewer:
  """
    Clase que genera una instancia de un objeto que contiene información de las Features (Dominio y Rangos) de una proteína según la base Uniprot pasada por su código pdb.
    Tiene métodos para gráficar a la proteína y sus features coloreandolos de forma particular.
  """
  def __init__(self, codigo_pdb, dominios=None, regiones=None, colores=None):
    self.codigo_pdb = codigo_pdb
    self.features = None #DataFrame con features
    self.dominios = dominios #Diccionario con dominios
    self.regiones = regiones #Diccionario con regiones
    self.colores = colores #Por si el usuario quiere "pintar" manualmente los dominios o regiones

    # Ejecución automática una vez inicializada la clase.
    self.busqueda_features_uniprot(self.codigo_pdb)
    if self.features is not None and not self.features.empty:
      self.obtener_dominios(self.features.iloc[0]['Domain [FT]'])
      self.obtener_regiones(self.features.iloc[0]['Region'])


  def busqueda_features_uniprot(self, pdb):
    """
      Realiza una búsqueda a la base UniProtKB según codigo PDB para obtener las Features de Dominios y Regiones de la proteína.
      Retorna un tabla de anotaciones.
    """
    consulta = f"(xref:pdb-{pdb})"
    parametros = {
      "query": consulta,
      "format": "tsv",
      "fields": "accession,id,gene_primary,ft_domain,ft_region"
    }
    busquedaUniprot = API_Uniprot_rest(parametros, URL_Uniprot)
    df_pdb = busquedaUniprot.data_de_paginacion_tsv(1)
    self.features = df_pdb
    return df_pdb

  def obtener_dominios(self, texto_dominio_uniprot):
    """
    Parsea una cadena con dominios y notas embebidas en una sola línea.
    """
    patron = r'DOMAIN (\d+)\.\.(\d+); /note="(.*?)"'
    matches = re.findall(patron, texto_dominio_uniprot)
    dominios = [{'inicio': int(start), 'final': int(end), 'note': note} for start, end, note in matches]
    self.dominios = dominios
    return dominios

  def obtener_regiones(self, texto_region_uniprot):
    """
    Parsea una cadena con regiones y notas embebidas en una sola línea.
    """
    patron = r'REGION (\d+)\.\.(\d+); /note="(.*?)"'
    matches = re.findall(patron, texto_region_uniprot)
    regiones = [{'inicio': int(start), 'final': int(end), 'note': note} for start, end, note in matches]
    self.regiones = regiones
    return regiones

  def nombre_color_masCercano(self, color_hex):
    """Devuelve el nombre del color CSS más cercano al color dado en formato #RRGGBB."""
    def hex_to_rgb(hex):
        hex = hex.lstrip('#')
        return tuple(int(hex[i:i+2], 16) for i in (0, 2 ,4))

    r1, g1, b1 = hex_to_rgb(color_hex)
    distancia_minima = float('inf')
    nombre_cercano = None

    for nombre, valor_hex in CSS4_COLORS.items():
        r2, g2, b2 = hex_to_rgb(valor_hex)
        dist = (r1 - r2)**2 + (g1 - g2)**2 + (b1 - b2)**2
        if dist < distancia_minima:
            distancia_minima = dist
            nombre_cercano = nombre
    return nombre_cercano

  def random_color(self):
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

  def mostrar_pdb_domains_regiones(self):
    # Genera la vista
    view = py3Dmol.view(query='pdb:' + self.codigo_pdb)
    view.setStyle({'cartoon': {'color': 'lightgrey'}})

    leyenda_dominios, leyenda_regiones = [], []
    # Dominios
    if self.dominios:
      for dominio in self.dominios:
        color = self.random_color()
        residuo = list(range(dominio['inicio'], dominio['final']+1))
        view.addStyle({'resi': residuo}, {'cartoon': {'color': color}})
        colorNombre = self.nombre_color_masCercano(color)
        leyenda_dominios.append((dominio['note'], dominio['inicio'], dominio['final'], color, colorNombre))
        
    # Regiones
    if self.regiones:
      for region in self.regiones:
        color = self.random_color()
        residuo = list(range(region['inicio'], region['final']+1))
        view.addStyle({'resi': residuo}, {'cartoon': {'color': color}})
        colorNombre = self.nombre_color_masCercano(color)
        leyenda_regiones.append((region['note'], region['inicio'], region['final'], color, colorNombre))

    view.zoomTo()
    html = view._make_html()
    match = re.search(r'<body>(.*?)</body>', html, flags=re.DOTALL)
    cuerpo = match.group(1) if match else html

    html_completo = self.generar_html_completo(
      self.codigo_pdb, cuerpo, leyenda_dominios, leyenda_regiones
    )
    
    # Guardar y abrir
    html_nombre = f"{self.codigo_pdb}_estructura_pdb.html"
    with open(html_nombre, 'w', encoding='utf-8') as f:
      f.write(html_completo)
    abrir_en_navegador(html_nombre)
    print(f"{html_nombre} guardado y abierto exitosamente.")
  
  def generar_html_completo(self, codigo_pdb, cuerpo_contenido, leyenda_dominios, leyenda_regiones):
    """
    Genera y retorna el HTML completo de la página.
    """
    def construir_seccion(titulo, items):
        html = f"<section class='legend-section'><h2>{titulo}</h2><ul>"
        for nota, inicio, fin, color, colorNombre in items:
            html += (
                f"<li><span class='color-box' style='background:{color};'></span>"
                f"<strong>{nota} ({inicio}\u2013{fin})</strong>: {colorNombre} ({color})</li>"
            )
        html += '</ul></section>'
        return html

    leyenda_html = construir_seccion('Dominios', leyenda_dominios)
    leyenda_html += construir_seccion('Regiones', leyenda_regiones)

    # Template completo con estilo más sofisticado
    html = f"""<!DOCTYPE html>
<html lang='es'>
<head>
  <meta charset='utf-8'/>
  <meta name='viewport' content='width=device-width, initial-scale=1.0'>
  <title>Estructura de {codigo_pdb}</title>
  <style>
    :root {{
      --bg-color: #f5f5f5;
      --card-bg: #ffffff;
      --primary: #333333;
      --accent: #1e88e5;
      --border-radius: 8px;
      --box-shadow: 0 4px 8px rgba(0,0,0,0.1);
    }}
    * {{ box-sizing: border-box; margin:0; padding:0; }}
    body {{ background: var(--bg-color); color: var(--primary); font-family: 'Segoe UI', Tahoma, sans-serif; }}
    .container {{ max-width: 1200px; margin: auto; padding: 20px; display: flex; flex-direction: column; align-items: center; }}
    h1 {{ margin-bottom: 16px; color: var(--accent); }}
    .card {{ background: var(--card-bg); border-radius: var(--border-radius); box-shadow: var(--box-shadow); width: 100%; margin-bottom: 24px; padding: 16px; }}
    .visualization {{ height: 600px; }}
    .legend {{ display: flex; flex-direction: column; gap: 16px; }}
    .legend-section h2 {{ font-size: 1.2em; margin-bottom: 8px; border-bottom: 2px solid var(--accent); padding-bottom: 4px; }}
    .legend-section ul {{ list-style: none; }}
    .legend-section li {{ display: flex; align-items: center; margin-bottom: 6px; }}
    .legend-section .color-box {{ width: 16px; height: 16px; border-radius: 4px; margin-right: 8px; border: 1px solid #ccc; }}
    @media (max-width: 768px) {{
      .visualization {{ height: 400px; }}
      .card {{ padding: 12px; }}
    }}
  </style>
</head>
<body>
  <div class='container'>
    <h1>Estructura de {codigo_pdb}</h1>
    <div class='card visualization'>
      {cuerpo_contenido}
    </div>
    <div class='card legend'>
      {leyenda_html}
    </div>
  </div>
</body>
</html>"""
    return html
  
def abrir_en_navegador(html_file):
    ruta_absoluta = os.path.abspath(html_file)
    url = f'file://{ruta_absoluta}'

    # Detectar navegadores disponibles
    navegadores_preferidos = ['vivaldi', 'firefox', 'brave', 'chrome', 'google-chrome']

    navegador_disponible = None
    for navegador in navegadores_preferidos:
        path = shutil.which(navegador)
        if path:
            navegador_disponible = path
            break

    if navegador_disponible:
        try:
            webbrowser.get(f'"{navegador_disponible}" %s').open_new_tab(url)
            print(f"Abriendo con: {navegador_disponible}")
        except webbrowser.Error:
            print("Error al intentar abrir con navegador específico. Usando navegador por defecto.")
            webbrowser.open_new_tab(url)
    else:
        print("No se detectó Chrome o Firefox. Abriendo con navegador por defecto.")
        webbrowser.open_new_tab(url)