import py3Dmol
import re
import pandas as pd
import requests
import random
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS

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

class PDB_3D_Viwer:
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
    """
      Muestra la estructura proteica y sus dominios segun anotaciones de Uniprot (Se requiere texto parseado de los rangos de sus dominios y nombres.)
    """
    view = py3Dmol.view(query='pdb:' + self.codigo_pdb)
    view.setStyle({'cartoon': {'color': 'lightgrey'}})
    leyenda_dominios = []
    leyenda_regiones = []

    if self.dominios:
      for dominio in self.dominios:
        color = self.random_color()
        rango_dominio = list(range(dominio['inicio'], dominio['final'] + 1))
        view.addStyle({'resi': rango_dominio}, {'cartoon': {'color': color}})
        leyenda_dominios.append([dominio['note'], dominio['inicio'], dominio['final'], color])

    if self.regiones:
      for region in self.regiones:
        color = self.random_color()
        rango_region = list(range(region['inicio'], region['final'] + 1))
        view.addStyle({'resi': rango_region}, {'cartoon': {'color': color}})
        leyenda_regiones.append([region['note'], region['inicio'], region['final'], color])


    view.zoomTo()
    view.show()

     # Imprimir leyenda
    print("Leyenda de dominios:")
    for note, inicio, final, color in leyenda_dominios:
      color_nombre = self.nombre_color_masCercano(color)
      print(f"🔹 {note} ({inicio}–{final}) → {color_nombre} ({color})")

    print("Leyenda de Regiones:")
    for note, inicio, final, color in leyenda_regiones:
      color_nombre = self.nombre_color_masCercano(color)
      print(f"🔸 {note} ({inicio}–{final}) → {color_nombre} ({color})")
