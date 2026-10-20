"""Tablero Streamlit - Trabajo Final Bioinformatica I.
Ejecutar:  pip install -r requirements.txt  &&  streamlit run app.py
"""
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

TITULO = "Trabajo Final Introduccion a la Bioinformatica"
st.set_page_config(page_title="Trabajo Final Introduccion a la Bioinformatica", page_icon="🧬", layout="centered")
st.markdown("<style>.block-container{max-width:1050px;padding-top:2rem}"
            "div[data-testid='stMetricValue']{color:#ff6b6b}</style>", unsafe_allow_html=True)

AD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")
ROJO, VIOLETA, GRIS, AMBAR, PLASMA = "#e5383b", "#9d6fd0", "#a89a9a", "#f4a261", "#f2d49b"  # arterial / venosa / hemo-hierro / plasma
DD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def img(nombre, caption=None):
    ruta = os.path.join(AD, nombre + ".jpg")
    if os.path.exists(ruta):
        st.image(ruta, caption=caption, width="stretch")
    else:
        st.warning(f"Falta la imagen assets/{nombre}.jpg")


def estilo(fig, h=380):
    fig.update_layout(template="plotly_dark", height=h, paper_bgcolor="rgba(0,0,0,0)",
                      plot_bgcolor="rgba(0,0,0,0)", margin=dict(l=10, r=10, t=40, b=10))
    return fig



import json
import urllib.request
import streamlit.components.v1 as components


def pdb_texto(pid):
    """Lee data/<ID>.pdb si existe; si no, lo baja de RCSB (y lo guarda ahí)."""
    for nombre in (f"{pid}.pdb", f"{pid.lower()}.pdb"):
        ruta = os.path.join(DD, nombre)
        if os.path.exists(ruta):
            return open(ruta, encoding="utf-8", errors="ignore").read()
    try:
        with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=25) as r:
            txt = r.read().decode("utf-8", "ignore")
        os.makedirs(DD, exist_ok=True)
        open(os.path.join(DD, f"{pid}.pdb"), "w", encoding="utf-8").write(txt)
        return txt
    except Exception:
        return None


_JS3D = r"""
const v = $3Dmol.createViewer("v", {backgroundColor: "#12080a"});
v.addModel(__PDB__, "pdb");
v.setStyle({chain: ["A","C"]}, {cartoon: {color: "#8d7f86"}});
v.setStyle({chain: ["B","D"]}, {cartoon: {color: "#d9a5a5"}});
if (__HEMO__) {
  v.addStyle({resn: "HEM"}, {stick: {color: "#f4a261", radius: 0.2}});
  v.addStyle({resn: "HEM", atom: "FE"}, {sphere: {color: "#ff9f1c", radius: 0.9}});
}
if (__CONTR__) {
  v.addStyle({chain: ["B","D"], resi: [42,63,67,92,98]}, {stick: {color: "#9d6fd0", radius: 0.25}});
  v.addStyle({chain: ["B","D"], resi: 6}, {sphere: {color: "#f2d49b"}});
}
v.addStyle({chain: ["B","D"], resi: 26}, {sphere: {color: "#ff2d2d"}});
v.addResLabels({chain: ["B","D"], resi: 26}, {fontColor: "white", backgroundColor: "#7a1018", backgroundOpacity: 0.8, fontSize: 13});
if (__DIST__) {
  const ca = v.selectedAtoms({chain: "B", resi: 26, atom: "CA"})[0];
  const fe = v.selectedAtoms({resn: "HEM", atom: "FE", chain: "B"})[0];
  if (ca && fe) {
    const d = Math.sqrt((ca.x-fe.x)**2 + (ca.y-fe.y)**2 + (ca.z-fe.z)**2);
    v.addLine({start: ca, end: fe, dashed: true, color: "white"});
    v.addLabel(d.toFixed(1) + " Å (Cα–Fe)", {position: {x:(ca.x+fe.x)/2, y:(ca.y+fe.y)/2, z:(ca.z+fe.z)/2},
      fontColor: "white", backgroundColor: "#000", backgroundOpacity: 0.6, fontSize: 12});
  }
}
if (__SURF__) v.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.45, color: "#c9b3b3"}, {chain: ["B"]});
v.zoomTo(); v.render();
if (__SPIN__) v.spin("y", 0.6);
"""


def visor3d(pid, hemo=True, contraste=False, dist=False, superficie=False, girar=False, alto=520):
    txt = pdb_texto(pid)
    if txt is None:
        st.caption(f"Visor 3D no disponible por ahora: no se pudo obtener {pid} (se descarga de RCSB al abrir el tablero con internet).")
        return
    lib = open(os.path.join(AD, "3Dmol-min.js"), encoding="utf-8").read()
    js = (_JS3D.replace("__PDB__", json.dumps(txt)).replace("__HEMO__", str(hemo).lower())
          .replace("__CONTR__", str(contraste).lower()).replace("__DIST__", str(dist).lower())
          .replace("__SURF__", str(superficie).lower()).replace("__SPIN__", str(girar).lower()))
    html = (f'<div id="v" style="width:100%;height:{alto}px;position:relative;border-radius:10px;overflow:hidden"></div>'
            f"<script>{lib}</script><script>{js}</script>")
    components.html(html, height=alto + 8)


# ---------------------------------------------------------------- datos del notebook
PANEL = pd.DataFrame([
    ("CD41/42", "HBB:c.126_129del", "Frameshift"), ("CD17", "HBB:c.52A>T", "Nonsense"),
    ("IVS-II-654", "HBB:c.316-197C>T", "Splicing"), ("-28", "HBB:c.-78A>G", "Promoter"),
    ("HbE", "HBB:c.79G>A", "Missense"), ("-50", "HBB:c.-100G>A", "Promoter"),
    ("Hb NewYork", "HBB:c.341T>A", "Missense"), ("CD71/72", "HBB:c.217dupA", "Frameshift"),
    ("IVS-I-1", "HBB:c.92+1G>T", "Splicing"), ("CD43", "HBB:c.130G>T", "Nonsense"),
    ("-29", "HBB:c.-79A>G", "Promoter"), ("CD27/28", "HBB:c.85dupC", "Frameshift"),
    ("IVS-II-5", "HBB:c.315+5G>C", "Splicing")], columns=["Mutación", "HGVS", "Consecuencia molecular"])

POS = pd.DataFrame([  # Sección 3.4.1 (Biopython sobre 2HHB, cadena B)
    ("β26 · HbE (Glu→Lys)", 26, 14.3, 6.1, 13.6, 19.0, "Objetivo", "Leve/moderada (β⁺-tal)"),
    ("β113 · Hb New York", 113, 16.4, 5.6, 13.8, 2.4, "Panel (leve)", "Leve"),
    ("β92 · His F8 proximal", 92, 2.1, 19.3, 5.6, 0.0, "Control (coordina el Fe)", "Control"),
    ("β63 · His E7 distal", 63, 4.1, 15.7, 12.6, 0.0, "Bolsillo del hemo", "Control"),
    ("β42 · Hb Hammersmith", 42, 6.7, 14.9, 7.2, 0.0, "Bolsillo del hemo", "Severa (inestable)"),
    ("β67 · Hb Bristol", 67, 4.1, 14.9, 14.2, 4.8, "Bolsillo del hemo", "Severa (inestable)"),
    ("β98 · Hb Köln", 98, 5.8, 16.9, 4.1, 0.0, "Bolsillo del hemo", "Severa (inestable)"),
    ("β6 · HbS (Glu→Val)", 6, 27.8, 12.1, 29.7, 109.4, "Superficie", "Severa (polimerización)")],
    columns=["Variante", "Residuo", "d Fe (Å)", "d α1β1 (Å)", "d α1β2 (Å)", "SASA (Å²)", "Contexto", "Severidad"])

# Hélice B de la β-globina madura (UniProt P68871, sin Met): residuos 20-32
SEQ_B = "VDEVGGEALGRLL"
# ---------------------------------------------------------------- barra lateral
secciones = [("s0", "Resumen"), ("s1", "1 · Mutación y hemoglobina"), ("s2", "2 · La paradoja de HbE"),
             ("s3", "3 · HbE vs. malaria y mapa"), ("s4", "4 · El panel de 13 variantes"),
             ("s5", "5 · Metodología"), ("s6", "6 · Conservación evolutiva"),
             ("s7", "7 · Mapeo estructural 3D"), ("s8", "8 · Cristal HbE vs. silvestre"),
             ("s9", "9 · Integración"), ("s10", "10 · Conclusiones"), ("s11", "Glosario y bibliografía")]
with st.sidebar:
    st.markdown("### 🧬 Índice")
    for a, t in secciones:
        st.markdown(f"[{t}](#{a})")
    st.caption("Bioinformática 2026")

# ---------------------------------------------------------------- encabezado
st.markdown("""<div class="hero"><svg viewBox="0 0 100 100" width="74" height="74"><defs>
<radialGradient id="g" cx="50%" cy="50%" r="50%"><stop offset="0%" stop-color="#7a1018"/><stop offset="45%" stop-color="#c81d25"/>
<stop offset="100%" stop-color="#ff5a5f"/></radialGradient></defs><ellipse cx="50" cy="50" rx="46" ry="40" fill="url(#g)"/>
<ellipse cx="50" cy="50" rx="17" ry="12" fill="#7a1018" opacity=".55"/></svg>
<div><div class="t">Hemoglobina · β-globina · HBB</div></div></div>""", unsafe_allow_html=True)
st.title(TITULO)
st.caption("Francisco Jose Rizzi · Nicolas Wald · Victor Orsingher · La paradoja de la Hemoglobina E: análisis estructural-mecanístico "
           "in silico de la variante β-globina más prevalente del Sudeste Asiático · a partir del panorama evolutivo de Deng et al. (2025)")

st.header("Resumen en 30 segundos", anchor="s0")
st.info("**La pregunta:** ¿por qué una mutación de tan alto impacto clínico y tan alta frecuencia (la Hemoglobina E) "
        "es, vista como proteína, casi benigna?  \n"
        "**La respuesta de este trabajo:** el cambio de aminoácido (Glu→Lys en la posición 26 de la β-globina) cae en un "
        "lugar de la molécula que tolera el cambio. El daño real de HbE no está en la forma de la proteína sino en el "
        "ARN mensajero (activa un sitio de *splicing* críptico). Que la proteína tolere el cambio explica que la "
        "variante pudiera hacerse muy frecuente bajo la presión de la malaria.")
c1, c2, c3, c4 = st.columns(4)
c1.metric("Distancia β26 → hierro del hemo", "14,3 Å", help="Mínima interatómica en 2HHB (Biopython)")
c2.metric("Desplazamiento T↔R de β26", "0,31 Å")
c3.metric("RMSD cristal HbE vs. silvestre", "0,525 Å", help="Cα de la cadena β, 3DUT vs 2HHB")
c4.metric("Glu en β26 (mamíferos)", "97,4 %", help="221 de 227 ortólogos de β-globina")

# ---------------------------------------------------------------- 1
st.header("1 · Para empezar: ¿qué es la hemoglobina y qué es una mutación?", anchor="s1")
st.markdown("""
La **hemoglobina** es la proteína de los glóbulos rojos que transporta oxígeno. Es un tetrámero: dos cadenas **α** y dos **β**,
y cada una guarda un **grupo hemo** con un átomo de hierro que es el que atrapa el O₂. La cadena β la fabrica el gen **HBB**.

Una **mutación missense** cambia una letra del ADN y, con ella, un aminoácido de la proteína. En la **Hemoglobina E (HbE)**
la letra **G→A** (HBB:c.79G>A) convierte el codón GAG en AAG: el **ácido glutámico (Glu, carga −1)** pasa a **lisina (Lys, carga +1)**
en la posición **β26**. Es una inversión de carga completa.
""")
col1, col2 = st.columns([1, 1.3])
with col1:
    forma = st.radio("Elegí la versión del gen", ["Silvestre (HbA)", "Hemoglobina E"], horizontal=True)
    es_e = forma.startswith("Hemo")
    st.metric("Codón 26 (proteína madura)", "AAG" if es_e else "GAG")
    st.metric("Aminoácido", "Lisina (Lys, K)" if es_e else "Glutámico (Glu, E)")
    st.metric("Carga a pH fisiológico", "+1" if es_e else "−1")
with col2:
    seq = list(SEQ_B)
    pos = [f"{20 + i}" for i in range(len(seq))]
    if es_e:
        seq[6] = "K"
    colores = [ROJO if p == "26" else (AMBAR if p == "30" else "#3a2026") for p in pos]
    fig = go.Figure(go.Bar(x=pos, y=[1] * len(seq), marker_color=colores, text=seq, textposition="inside",
                           textfont=dict(size=20, color="white"), hoverinfo="skip"))
    fig.update_yaxes(visible=False)
    fig.update_xaxes(title="Posición en la β-globina madura (hélice B, residuos 20–32)", type="category")
    fig.update_layout(title="Entorno de β26 (rojo). En ámbar, Arg30, su compañera de puente salino en el estado R")
    st.plotly_chart(estilo(fig, 300), width="stretch")
st.caption("Nota de numeración: las estructuras cristalográficas usan la proteína madura (sin la Met inicial), así β26 = residuo 26. "
           "AlphaFold y las secuencias completas incluyen la Met, y entonces la misma posición es el residuo 27.")

# ---------------------------------------------------------------- 2
st.header("2 · La paradoja de la Hemoglobina E", anchor="s2")
st.markdown("""
Deng et al. (2025) reconstruyeron el origen y la historia de las mutaciones de **β-talasemia** en el sur de China, la enfermedad
monogénica más prevalente del Sudeste Asiático. Trece variantes explican ≈95 % de los casos; la mayoría rompe la **expresión** del gen
y solo dos son *missense*. Una de ellas, **HbE**, es doblemente peculiar:

- es una **hemoglobina anómala** (cambia un aminoácido), **y**
- es un **alelo β⁺-talasémico**: la sustitución activa un **sitio críptico de *splicing*** en el exón 1 (Orkin et al., 1982), de modo que
  se produce menos β-globina normal.

Participa en cerca del 50 % de la talasemia severa a nivel mundial. Deng et al. la tratan como un problema poblacional y evolutivo, **sin
analizar el mecanismo a nivel de la proteína**. Ese es el hueco que cubre este trabajo.
""")
st.success("**Hipótesis:** β26 ocupa una posición superficial (B8), alejada del hemo y de la interfaz alostérica α1β2. "
           "Esperamos hallar una posición estructuralmente permisiva y, en contraste, que las hemoglobinas inestables severas "
           "mapeen en posiciones enterradas, de contacto con el hemo o de interfaz.")
st.markdown("**Resultados esperados:** (i) clasificación mecanística del panel; (ii) mapa 3D y perfil de conservación de la posición; "
            "(iii) comparación silvestre vs. mutante; (iv) contraste con variantes estructurales severas.")

# ---------------------------------------------------------------- 3
st.header("3 · HbE y malaria: por qué una mutación así llegó a ser tan frecuente", anchor="s3")
st.markdown("""
En 1949 Haldane propuso que ciertas variantes de hemoglobina se mantienen en la población porque **los portadores (heterocigotas)
resisten mejor la malaria**. El caso clásico es HbS en África. En el Sudeste Asiático el equivalente es HbE: está concentrada donde la
malaria por *Plasmodium falciparum* fue endémica, con su máximo en la zona donde se juntan **Tailandia, Laos y Camboya** y una frecuencia que
decae en círculos concéntricos hacia afuera (Nagel et al., 1981).
""")

st.subheader("3.1 · ¿Qué tan frecuente es HbE?")
freq = pd.DataFrame({"Población": ["So (Khammuan, Laos)", "Alak/Ngeh (Sekong, Laos)", "Oy (Attapeu, Laos)",
                                   "NE de Tailandia (muestras grandes, ≈)"], "q": [0.426, 0.433, 0.253, 0.30]})
vista = st.radio("Qué mostrar", ["Frecuencia del alelo HbE (q)", "Proporción esperada de cada genotipo (Hardy–Weinberg)"], horizontal=True)
if vista.startswith("Frecuencia"):
    f = px.bar(freq, x="q", y="Población", orientation="h", text="q", color_discrete_sequence=[ROJO])
    f.update_traces(texttemplate="%{text:.3f}")
    f.update_xaxes(title="Frecuencia del alelo (0–1)", range=[0, 0.5])
else:
    d = freq.assign(Silvestres=(1 - freq.q) ** 2, **{"Portadores (AE)": 2 * freq.q * (1 - freq.q)}, Homocigotas=freq.q ** 2)
    d = d.melt(id_vars="Población", value_vars=["Silvestres", "Portadores (AE)", "Homocigotas"], var_name="Genotipo", value_name="Proporción")
    f = px.bar(d, x="Proporción", y="Población", color="Genotipo", orientation="h", barmode="stack",
               color_discrete_map={"Silvestres": GRIS, "Portadores (AE)": VIOLETA, "Homocigotas": ROJO})
st.plotly_chart(estilo(f, 330), width="stretch")
st.caption("Fuente: Flatz et al., *Hemoglobin* 2004;28:197-204 (≈0,3 en el NE de Tailandia; 0,426 / 0,433 / 0,253 en tres poblaciones del sur de Laos). "
           "Con q≈0,43 casi la mitad de las personas es portadora. En el NE de la India se reportan portadores del 27–51 % en Assam/Meghalaya "
           "y 5–35 % en Bengala Occidental (OMIM).")

st.subheader("3.2 · Simulador: ventaja del heterocigota (modelo ilustrativo)")
st.markdown("Si los portadores (AE) son los más aptos, la frecuencia del alelo no desaparece ni llega al 100 %: se estabiliza en "
            "**q\\* = s / (s + t)**. Probá valores.")
s1, s2, s3, s4 = st.columns(4)
s = s1.slider("s · desventaja de los no portadores (malaria)", 0.0, 0.4, 0.10, 0.01)
t = s2.slider("t · costo del alelo en homocigotas / dobles heterocigotas", 0.01, 0.6, 0.20, 0.01)
q0 = s3.slider("Frecuencia inicial q₀", 0.001, 0.5, 0.01, 0.001, format="%.3f")
G = s4.slider("Generaciones", 50, 600, 300, 10)
q, tray = q0, [q0]
for _ in range(G):
    p = 1 - q
    wm = p * p * (1 - s) + 2 * p * q + q * q * (1 - t)
    q = (p * q + q * q * (1 - t)) / wm
    tray.append(q)
qeq = s / (s + t) if (s + t) > 0 else 0
f = go.Figure(go.Scatter(y=tray, mode="lines", line=dict(color=ROJO, width=3), name="q(t)"))
f.add_hline(y=qeq, line_dash="dash", line_color=VIOLETA, annotation_text=f"equilibrio q* = {qeq:.2f}")
f.update_xaxes(title="Generaciones"); f.update_yaxes(title="Frecuencia del alelo HbE", range=[0, 1])
st.plotly_chart(estilo(f, 340), width="stretch")
st.caption("Modelo de libro de texto (aptitudes: AA = 1−s, AE = 1, EE = 1−t), no un ajuste a datos reales. Con s=0,10 y t=0,20 el equilibrio (≈0,33) "
           "es del orden de lo observado.")

st.subheader("3.3 · ¿Qué evidencia hay de que la malaria seleccionó HbE?")
e1, e2 = st.columns(2)
with e1:
    st.markdown("**A favor**\n- *P. falciparum* se multiplica más lento en glóbulos rojos HbE/HbE (Nagel et al., 1981).\n"
                "- Los portadores tendrían menos parasitemias altas, coherente con la hipótesis de Haldane (Chotivanich et al., 2002, vía OMIM).\n"
                "- Desequilibrio de ligamiento extendido alrededor de HbE, huella típica de selección positiva (PMC1182083).\n"
                "- Máxima frecuencia justo donde la malaria fue hiperendémica.")
with e2:
    st.markdown("**Matices**\n- Un metaanálisis de dos estudios con HbE no encontró protección clara; otros trabajos sí (metaanálisis de hemoglobinopatías discutido en PMC4726726).\n"
                "- Una revisión concluye que hace falta más investigación antes de afirmar selección positiva (PMC6925914).\n"
                "- La protección parece frenar las parasitemias pesadas más que la infección leve.")
st.warning("La relación con la malaria **no se demuestra con las herramientas de este trabajo**: se apoya en la bibliografía (Thein, 2013). "
           "Lo que sí aporta el análisis estructural es la condición necesaria: que la proteína tolere el cambio.")

st.subheader("3.4 · Un fósil genético: la malaria de hoy ya no es la de entonces")
m1, m2, m3 = st.columns(3)
m1.metric("Casos de malaria, Gran Mekong 2000→2020", "−97 %")
m2.metric("Muertes estimadas por malaria (6 países)", "6000 → 15", help="2000 vs 2020, WHO")
m3.metric("Casos 2012→2022 (GMS)", "−77 %")
st.caption("Fuente: WHO (Día Mundial del Paludismo 2021) y Manzoni et al., *Malaria Journal* 2024. HbE persiste como huella de la presión pasada: "
           "los alelos tardan generaciones en cambiar de frecuencia. Cuidado: esto compara presiones en el tiempo, no es una correlación por individuo.")


st.subheader("3.5 · Mapa: ¿dónde están HbE y Hb New York, y dónde estuvo la malaria?")
st.markdown("""
Si la malaria sostenida en el tiempo empujó la frecuencia de HbE, esperaríamos que **donde HbE es más frecuente haya habido (histórica) malaria por
*P. falciparum*** intensa. Hb New York (β113 **Val→Glu**, ligeramente inestable, hallada en muchas familias chinas) sirve de contraste:
es una variante estructural **frecuente en el sur de China pero con prevalencia ~100 veces menor**.
""")
st.caption("Ojo: según las fuentes consultadas el máximo de HbE está en el **noreste de Tailandia, sur de Laos y Camboya**, no en el sur de Tailandia.")

HBE = pd.DataFrame([  # nombre, lat, lon, tipo, valor, fuente
    ("So (Khammuan, Laos)", 17.4, 105.2, "q", 0.426, "Flatz et al. 2004"),
    ("Alak/Ngeh (Sekong, Laos)", 15.5, 106.9, "q", 0.433, "Flatz et al. 2004"),
    ("Oy (Attapeu, Laos)", 14.8, 107.0, "q", 0.253, "Flatz et al. 2004"),
    ("NE de Tailandia (muestras grandes)", 15.5, 103.0, "q", 0.30, "Flatz et al. 2004 (≈)"),
    ("Assam/Meghalaya (NE de la India)", 25.8, 92.0, "port", 39.0, "OMIM/ClinVar (27–51 %, punto medio)"),
    ("Bengala Occidental (India)", 23.0, 87.9, "port", 20.0, "OMIM/ClinVar (5–35 %, punto medio)"),
    ("Chittagong Hill Tracts (Bangladesh, tribal)", 22.4, 92.2, "port", 41.7, "IthaMaps"),
    ("Bangladesh (país)", 23.7, 90.4, "port", 6.1, "IthaMaps")],
    columns=["Lugar", "lat", "lon", "tipo", "valor", "Fuente"])
HBE["Portadores HbE (%)"] = HBE.apply(lambda r: round(200 * r.valor * (1 - r.valor), 1) if r.tipo == "q" else r.valor, axis=1)
HNY = pd.DataFrame([
    ("Guangxi (tamizaje, n=30 691)", 22.8, 108.3, 0.12, "Chin J Hematol 2013"),
    ("Laibin, centro de Guangxi", 23.7, 109.2, 0.21, "Derivado: 0,68 % portadores de variantes × 30,3 % Hb New York"),
    ("Norte de Guangdong (Shaoguan)", 24.8, 113.6, 0.06, "Estudio de variantes en el norte de Guangdong (PMC10545308)")],
    columns=["Lugar", "lat", "lon", "Prevalencia (%)", "Fuente"])

LAT0, LAT1, LON0, LON1 = 5, 32, 84, 124


@st.cache_data
def cargar_malaria():
    """Capa de malaria incluida en data/malaria_pfpr.csv (columnas lat, lon, pfpr_2000, pfpr_2025; PfPR2-10 en %)."""
    ruta = os.path.join(DD, "malaria_pfpr.csv")
    return pd.read_csv(ruta) if os.path.exists(ruta) else None


def pfpr_sitio(df, lat, lon, radio=1.0):
    m = ((df.lat - lat) ** 2 + (df.lon - lon) ** 2) ** 0.5 <= radio
    return df.loc[m, "pfpr_2000"].mean(), df.loc[m, "pfpr_2025"].mean()


malaria = cargar_malaria()
capa = st.radio("Variantes a mostrar", ["HbE", "Hb New York", "Ambas"], horizontal=True, index=2)
anio = st.radio("Año de la capa de malaria", [2000, 2025], horizontal=True,
                format_func=lambda a: f"{a} ({'el más antiguo disponible' if a == 2000 else 'actual'})")
fm = go.Figure()
if malaria is not None:
    col = f"pfpr_{anio}"
    mm = malaria[malaria[col] >= 0.3]
    fm.add_trace(go.Scattergeo(lat=mm.lat, lon=mm.lon, mode="markers", name=f"Malaria (PfPR {anio})", opacity=0.8,
                               marker=dict(symbol="square", size=4, color=mm[col], cmin=0, cmax=30,
                                           colorscale=[[0, "#16302f"], [0.4, "#2ec4b6"], [1, "#e8fff9"]],
                                           colorbar=dict(title="PfPR (%)", len=0.5, x=0.02, thickness=12)),
                               hovertemplate="PfPR: %{marker.color:.1f} %<extra></extra>"))
if capa in ("HbE", "Ambas"):
    fm.add_trace(go.Scattergeo(lat=HBE.lat, lon=HBE.lon, mode="markers", name="HbE (% portadores)",
                               marker=dict(color=ROJO, size=8 + HBE["Portadores HbE (%)"] * 0.55, line=dict(color="white", width=1)),
                               text=HBE.Lugar + " · " + HBE.Fuente, customdata=HBE["Portadores HbE (%)"],
                               hovertemplate="%{text}<br>Portadores: %{customdata} %<extra></extra>"))
if capa in ("Hb New York", "Ambas"):
    fm.add_trace(go.Scattergeo(lat=HNY.lat, lon=HNY.lon, mode="markers", name="Hb New York (% población)",
                               marker=dict(color=VIOLETA, size=8 + HNY["Prevalencia (%)"] * 25, symbol="diamond", line=dict(color="white", width=1)),
                               text=HNY.Lugar + " · " + HNY.Fuente, customdata=HNY["Prevalencia (%)"],
                               hovertemplate="%{text}<br>Prevalencia: %{customdata} %<extra></extra>"))
fm.update_geos(projection_type="natural earth", lonaxis_range=[LON0, LON1], lataxis_range=[LAT0, LAT1], showland=True, landcolor="#1f1114",
               showocean=True, oceancolor="#12080a", showcountries=True, countrycolor="#5a3038", coastlinecolor="#7a4650", bgcolor="rgba(0,0,0,0)")
fm.update_layout(legend=dict(orientation="h", y=-0.02))
st.plotly_chart(estilo(fm, 520), width="stretch")
st.caption("Escalas distintas: HbE = % de portadores (para q se usó 2q(1−q)); Hb New York = % de la población tamizada. Los estudios usan métodos y poblaciones diferentes, "
           "así que sirve para comparar **zonas**, no valores exactos.")
if malaria is not None:
    filas = []
    for nombre, la, lo in list(HBE[["Lugar", "lat", "lon"]].itertuples(index=False)) + list(HNY[["Lugar", "lat", "lon"]].itertuples(index=False)):
        a0, a1 = pfpr_sitio(malaria, la, lo)
        filas.append((nombre, round(a0, 2), round(a1, 2)))
    filas.append(("Promedio de toda la región mostrada", round(malaria.pfpr_2000.mean(), 2), round(malaria.pfpr_2025.mean(), 3)))
    st.markdown("**Malaria alrededor de cada sitio** (PfPR2-10 medio en un radio de ~110 km, %)")
    st.dataframe(pd.DataFrame(filas, columns=["Sitio", "PfPR 2000 (%)", "PfPR 2025 (%)"]), hide_index=True, width="stretch")
    st.info("**Lectura honesta del mapa.** En 2000 la malaria modelada es claramente mayor que el promedio regional (0,67 %) en el sur de Laos "
            "(donde HbE llega a q≈0,43) y se mantiene en torno a Chittagong (Bangladesh). En cambio, en el noreste de Tailandia y en el sur de China "
            "la parasitemia ya era casi nula en 2000. Eso no refuta la hipótesis: estas capas son una foto **posterior** a décadas de control, "
            "y la selección que dejó HbE actuó durante siglos. Es decir, el mapa actual no reconstruye la presión pasada; "
            "para eso haría falta cartografía histórica de endemicidad. Por eso 3.4 habla de HbE como un fósil genético.")
with st.expander("Datos del mapa y fuentes"):
    st.dataframe(HBE[["Lugar", "Portadores HbE (%)", "Fuente"]], hide_index=True, width="stretch")
    st.dataframe(HNY[["Lugar", "Prevalencia (%)", "Fuente"]], hide_index=True, width="stretch")
st.warning("**Qué se puede y qué no se puede concluir.** La coincidencia geográfica entre HbE y malaria es coherente con la hipótesis de selección, "
           "pero es una asociación entre zonas, no una prueba individual (falacia ecológica), y la malaria actual es mucho menor que la histórica. "
           "Tampoco la bioinformática de este trabajo demuestra la selección: aporta que la proteína **tolera** el cambio. "
           "La evidencia directa viene de estudios de laboratorio y genética de poblaciones (sección 3.3), que a la fecha no son unánimes.")

# ---------------------------------------------------------------- 4
st.header("4 · Triage mecanístico: 11 de 13 variantes ni siquiera llegan a la proteína", anchor="s4")
st.markdown("Antes de mirar estructuras, se clasificaron las trece variantes del panel de Deng et al. por su consecuencia molecular "
            "(`SupplementaryData1.xlsx`). Solo las **missense** producen una proteína completa con *una* sustitución que se pueda analizar en 3D.")
cuenta = PANEL["Consecuencia molecular"].value_counts().reset_index()
cuenta.columns = ["Consecuencia", "Variantes"]
cuenta["Nivel"] = cuenta["Consecuencia"].map(lambda x: "Estructura (proteína)" if x == "Missense" else "Expresión del gen")
f = px.bar(cuenta, x="Consecuencia", y="Variantes", color="Nivel", text="Variantes",
           color_discrete_map={"Expresión del gen": VIOLETA, "Estructura (proteína)": ROJO})
st.plotly_chart(estilo(f, 330), width="stretch")
sel = st.multiselect("Filtrar la tabla por consecuencia", sorted(PANEL["Consecuencia molecular"].unique()),
                     default=sorted(PANEL["Consecuencia molecular"].unique()))
st.dataframe(PANEL[PANEL["Consecuencia molecular"].isin(sel)], hide_index=True, width="stretch")
st.markdown("Se analizó estructuralmente solo **HbE** (alta frecuencia poblacional y doble carácter: hemoglobinopatía + alelo β⁺-talasémico); "
            "Hb New York entra como variante de contraste.")

# ---------------------------------------------------------------- 5
st.header("5 · Metodología: cómo se hizo", anchor="s5")
st.graphviz_chart("""digraph { rankdir=LR; bgcolor="transparent"; node [shape=box, style="rounded,filled", fillcolor="#231014", color="#e5383b", fontcolor="#f5ebe8", fontsize=11];
edge [color="#a89a9a"];
A [label="Panel Deng 2025\\n(13 variantes)"]; B [label="Triage\\nmecanístico"]; C [label="PSI-BLAST\\n(ortólogos β)"]; D [label="MAFFT +\\nconservación"];
E [label="Mapeo en 2HHB / 1HHO\\n(distancias, SASA, VMD)"]; F [label="AlphaFold\\nWT vs E26K"]; G [label="Cristal 3DUT vs 2HHB\\n(RMSD)"]; H [label="Integración"];
A->B; B->C; C->D; D->E; E->F; F->G; G->H; }""")
datos = pd.DataFrame([
    ("Secuencia β-globina humana", "UniProt P68871", "Secuencia de referencia (147 aa con Met)"),
    ("2HHB", "PDB", "Hb humana desoxigenada, estado T (silvestre)"),
    ("1HHO", "PDB", "Hb humana oxigenada, estado R (silvestre)"),
    ("3DUT", "PDB", "Hb E desoxigenada (Glu26Lys), estado T, tetrámero α₂β₂ completo con hemo"),
    ("Panel de variantes", "Deng et al. 2025 (suplementarios) + HbVar", "Catálogo de variantes"),
    ("Modelos AlphaFold", "WT y E26K", "Predicción de confianza y de cambio estructural")], columns=["Dato", "Origen", "Uso"])
st.dataframe(datos, hide_index=True, width="stretch")
st.caption("Herramientas: PSI-BLAST · MAFFT (L-INS-i) · pymsaviz · Biopython (PDBParser, ShrakeRupley, Superimposer) · VMD. "
           "Se compara 2HHB con 3DUT (ambas estado T) para aislar la mutación del reordenamiento alostérico T↔R.")

# ---------------------------------------------------------------- 6
st.header("6 · Conservación evolutiva de β26", anchor="s6")
st.markdown("Se buscaron ortólogos de β-globina en mamíferos con **PSI-BLAST** (RefSeq, 5 rondas de 500 secuencias, excluyendo α/γ/ε y *beta-like*, "
            "porque familias más amplias —neuroglobinas, no mamíferos— diluirían la señal). Se filtraron secuencias *predicted* y de baja calidad "
            "y se alinearon con **MAFFT**.")
with st.expander("Ver la búsqueda en PSI-BLAST"):
    a, b = st.columns(2)
    with a: img("psiblast_form", "Parámetros de la búsqueda")
    with b: img("psiblast_resultados", "Resultados (ronda final)")
cons = pd.DataFrame({"Residuo": ["E (Glu)", "H", "Q", "A", "R"], "Secuencias": [221, 2, 2, 1, 1]})
cons["%"] = (cons["Secuencias"] / 227 * 100).round(1)
f = px.bar(cons, x="Residuo", y="Secuencias", text="%", color="Residuo",
           color_discrete_sequence=[ROJO, GRIS, GRIS, GRIS, GRIS], log_y=True)
f.update_traces(texttemplate="%{text}%"); f.update_layout(showlegend=False, title="Aminoácido en la posición 26 (n = 227 β-globinas de mamíferos, eje log)")
st.plotly_chart(estilo(f, 340), width="stretch")
st.markdown("La posición es **esencialmente invariante: ácido glutámico en ≈97 %** de los ortólogos y **ninguna lisina** en ninguna especie. "
            "Las pocas excepciones (H, Q, A, R) coinciden con secuencias cuya longitud difiere de 147 aa, lo que desplaza el conteo "
            "(lectura a confirmar). Es decir: **la Lys de HbE no aparece de forma natural**.")
st.info("**Aparente contradicción central:** una posición tan conservada pero cuya sustitución es estructuralmente casi neutra. "
        "Se resuelve más abajo: la conservación no reflejaría una restricción de la proteína madura sino de otro nivel funcional.")
with st.expander("Ver el alineamiento múltiple (columna 27 = β26)"):
    img("msa", "Alineamiento MAFFT. La columna 27 corresponde a β26 porque el alineamiento incluye la Met inicial.")

# ---------------------------------------------------------------- 7
st.header("7 · Dónde está β26 en la estructura 3D", anchor="s7")
tabs = st.tabs(["Distancias y contexto", "Vista VMD", "Transición T↔R", "Puentes salinos", "Factor B y SASA", "AlphaFold"])
with tabs[0]:
    st.markdown("Para cada residuo de la cadena β de **2HHB** se midió la distancia mínima al hierro del hemo y a las interfaces α1β1 (empaquetamiento) "
                "y α1β2 (alostérica), más la accesibilidad al solvente (SASA, Shrake–Rupley).")
    st.dataframe(POS.drop(columns=["Residuo"]), hide_index=True, width="stretch")
    eje = st.selectbox("Distancia a mostrar", ["d Fe (Å)", "d α1β1 (Å)", "d α1β2 (Å)", "SASA (Å²)"])
    f = px.bar(POS.sort_values(eje), x=eje, y="Variante", orientation="h", color="Contexto", text=eje)
    st.plotly_chart(estilo(f, 420), width="stretch")
    f = px.scatter(POS, x="d Fe (Å)", y="SASA (Å²)", text="Variante", color="Severidad", size=[16] * len(POS))
    f.update_traces(textposition="top center")
    f.update_layout(title="Dos mapas en uno: lejos del hemo y expuesto ≠ inofensivo (compará HbE y HbS)")
    st.plotly_chart(estilo(f, 440), width="stretch")
    st.markdown("- β26 está a **14,3 Å del hemo** y a **13,6 Å de la interfaz α1β2**: lejos del sitio del O₂ y del mecanismo T↔R. "
                "A 6,1 Å de la interfaz α1β1 (sin contacto directo, que exigiría <4 Å).\n"
                "- Las variantes inestables severas (β42, β67, β98) y la His distal β63 están a **4–7 Å del Fe**, dentro del bolsillo del hemo. "
                "β92 (His proximal) a 2,1 Å **valida** el cálculo.\n"
                "- **HbS (β6)** es la posición más superficial (27,8 Å del Fe, SASA 109 Å²) y aun así es severa: su daño es la **polimerización "
                "intermolecular**, no un defecto local.")
with tabs[1]:
    st.markdown("Visor 3D interactivo (arrastrá para rotar, rueda para acercar). β26 de ambas cadenas β en **rojo**; el hemo en ámbar con el hierro "
                "como esfera. Las estructuras se bajan de RCSB (o se leen de `data/`).")
    ESTR = {"2HHB · Hb humana desoxigenada (estado T)": "2HHB", "1HHO · Hb humana oxigenada (estado R)": "1HHO",
            "3DUT · Hemoglobina E desoxigenada (Lys26)": "3DUT"}
    v1, v2 = st.columns([1.3, 1])
    elegida = v1.selectbox("Estructura", list(ESTR))
    o = v2.columns(2)
    hemo = o[0].checkbox("Hemos", True); contr = o[1].checkbox("Posiciones de contraste", False)
    o = v2.columns(3)
    dist = o[0].checkbox("Cα–Fe", False); sup = o[1].checkbox("Superficie β", False); gira = o[2].checkbox("Girar", False)
    visor3d(ESTR[elegida], hemo, contr, dist, sup, gira)
    st.caption("Posiciones de contraste: β42, β63, β67, β92, β98 en violeta (bolsillo del hemo) y β6 (HbS) en amarillo. "
               "La distancia mostrada es Cα–Fe, por eso difiere de la mínima interatómica de la tabla (14,3 Å).")
    with st.expander("Capturas originales de VMD del informe"):
        a, b, c = st.columns(3)
        with a: img("vmd_1")
        with b: img("vmd_2")
        with c: img("vmd_3")
with tabs[2]:
    st.markdown("Se alinearon 2HHB (desoxi, T) y 1HHO (oxi, R) por el esqueleto. β26 se desplaza solo **0,31 Å** entre estados: **no participa** del "
                "reordenamiento alostérico.")
    img("desplazamiento_T_R", "Desoxi en gris oscuro (β26 cian) · oxi en gris claro (β26 magenta)")
    d1, d2 = st.columns(2)
    with d1: img("dist_hemo_desoxi", "Estado T: 13,23 Å al hierro")
    with d2: img("dist_hemo_oxi", "Estado R: 14,13 Å al hierro")
    st.caption("Las mediciones en VMD (centro a centro) coinciden en orden de magnitud con la distancia mínima de Biopython (14,3 Å); difiere el criterio de medición.")
with tabs[3]:
    st.markdown("Se buscaron puentes salinos y de hidrógeno de la cadena β con VMD (corte O–N de 3,2 Å). **No hay enlaces de hidrógeno en β26** en ningún estado.")
    sb = pd.DataFrame({"Estado": ["Desoxi (2HHB, T)", "Oxi (1HHO, R)"], "Puentes salinos en cadena β": [6, 2],
                       "¿Involucra β26?": ["No", "Sí: Glu26–Arg30"]})
    st.dataframe(sb, hide_index=True, width="stretch")
    st.markdown("En el estado R, Glu26 forma un puente salino con Arg30 (misma hélice). Como E→K invierte la carga, **es la única interacción local "
                "que la sustitución podría perturbar**.")
with tabs[4]:
    st.markdown("**Factor B** (incertidumbre posicional en el cristal): β26 tiene un valor bajo, es decir, una posición bien determinada: "
                "la ubicación superficial es confiable y no un artefacto de una región desordenada.")
    img("factor_b", "Coloreado por factor B; β26 marcado con puntos blancos")
    st.markdown("**Accesibilidad al solvente (SASA)** de β26:")
    k1, k2, k3 = st.columns(3)
    k1.metric("Desoxi · VMD", "24,6 Å²"); k2.metric("Oxi · VMD", "37,8 Å²"); k3.metric("Desoxi · Biopython", "19,0 Å²")
    st.caption("Sonda de 1,4 Å. Mayor exposición en el estado oxigenado. Para dimensionar: β6 (HbS) tiene 109 Å² en el cálculo de Biopython, "
               "así que β26 está expuesta pero de forma moderada, no tanto como la posición más superficial del panel.")
    a, b = st.columns(2)
    with a: img("sasa_desoxi", "Superficie, estado desoxi")
    with b: img("sasa_oxi", "Superficie, estado oxi")
with tabs[5]:
    st.markdown("AlphaFold tiende a predecir el mismo plegamiento ante mutaciones puntuales, pero una caída de confianza puede delatar un efecto. "
                "Se compararon la silvestre y la mutante E26K (en AlphaFold la posición es el residuo 27).")
    a, b = st.columns(2)
    with a: img("af_silvestre", "Silvestre (β26 resaltado)")
    with b: img("af_mutante", "Mutante E26K")
    st.markdown("Ambos modelos tienen confianza muy alta y casi igual en el entorno de la mutación. Salvo en los extremos, el alineamiento no muestra "
                "diferencias mayores a 1 Å.")
    img("af_alineamiento", "Alineamiento estructural coloreado por RMSD")
    a, b = st.columns(2)
    with a: img("af_qres", "Qres por residuo")
    with b: img("af_rmsd", "RMSD por residuo")

# ---------------------------------------------------------------- 8
st.header("8 · La prueba más directa: el cristal real de HbE (3DUT vs. 2HHB)", anchor="s8")
st.markdown("AlphaFold es una predicción; **3DUT es la estructura experimental** de HbE. Ambas están en estado T, así que la diferencia se debe a la "
            "mutación y no a la transición alostérica. Se superpusieron las cadenas β por sus carbonos α (algoritmo de Kabsch). "
            "La identidad se verificó: posición 26 = Glu en 2HHB y **Lys en 3DUT**.")
rm = pd.DataFrame({"Comparación": ["Ruido interno: β1 vs β2 silvestres (2HHB)", "Global: β HbE vs β silvestre", "Local: hélice B (res. 20–32)"],
                   "RMSD Cα (Å)": [0.251, 0.525, 0.591]})
f = px.bar(rm, x="RMSD Cα (Å)", y="Comparación", orientation="h", text="RMSD Cα (Å)", color_discrete_sequence=[ROJO])
f.update_traces(texttemplate="%{text:.3f} Å"); f.add_vline(x=1, line_dash="dash", line_color=AMBAR, annotation_text="1 Å")
f.update_xaxes(range=[0, 1.2]); f.update_yaxes(autorange="reversed")
st.plotly_chart(estilo(f, 300), width="stretch")
st.success("Los tres valores están muy por debajo de 1 Å. La diferencia es del orden del doble del ruido interno y cae en el rango habitual entre "
           "dos cristales independientes de la misma proteína; el entorno de β26 no se distorsiona más que el promedio de la cadena. "
           "**La estructura experimental de HbE es prácticamente indistinguible de la silvestre.**")

# ---------------------------------------------------------------- 9
st.header("9 · Integración: dos posiciones superficiales, dos capas de daño", anchor="s9")
integ = pd.DataFrame([
    ("β26", "HbE (Glu→Lys)", "Superficie · 14,3 Å al hemo · 13,6 Å a α1β2", "Leve/moderada (β⁺-tal)", "Splicing críptico (capa de ARN)"),
    ("β6", "HbS (Glu→Val)", "Superficie · 27,8 Å al hemo", "Severa", "Polimerización intermolecular"),
    ("β42", "Hb Hammersmith", "Bolsillo del hemo · 6,7 Å al Fe", "Severa (inestable)", "Pérdida del hemo"),
    ("β67", "Hb Bristol", "Bolsillo del hemo · 4,1 Å al Fe", "Severa (inestable)", "Bolsillo del hemo alterado"),
    ("β98", "Hb Köln", "Enterrada/bolsillo · 5,8 Å al Fe", "Severa (inestable)", "Desestabilización del hemo")],
    columns=["Posición", "Variante", "Contexto en 2HHB", "Severidad", "Mecanismo del daño"])
st.dataframe(integ, hide_index=True, width="stretch")
st.markdown("En HbS una posición superficial produce enfermedad severa por polimerización. En HbE la permisividad de β26 (lejos del hemo y de α1β2, "
            "expuesta, casi sin movimiento en T↔R y con un único puente salino) es compatible con una función conservada, y el daño real está en el ARN. "
            "Esa tolerancia hace viable la variante y, por lo tanto, **seleccionable por la malaria** (Thein, 2013).")
b1, b2 = st.columns(2)
b1.markdown("**Con datos propios (medidos aquí)**\n- conservación de β26\n- distancias al hemo y a interfaces\n- SASA y factor B\n"
            "- desplazamiento T↔R\n- puentes salinos\n- RMSD de 3DUT vs 2HHB y de AlphaFold")
b2.markdown("**Apoyado en la bibliografía**\n- splicing críptico en el exón 1 (Orkin et al., 1982)\n- selección por malaria (Thein, 2013; sección 3 de este tablero)")

# ---------------------------------------------------------------- 10
st.header("10 · Conclusiones y proyección", anchor="s10")
st.markdown("""
1. La sustitución **Glu→Lys en β26 invierte la carga pero casi no altera la estructura**: la lisina queda en superficie, lejos del hemo y de la interfaz alostérica.
2. Modelos de AlphaFold y, sobre todo, el **cristal experimental 3DUT** (RMSD 0,525 Å) lo confirman.
3. La β26 está **muy conservada** en mamíferos y, sin embargo, tolera el cambio: la restricción evolutiva no vendría de la estructura de la proteína madura.
4. El impacto clínico de HbE **no se explica por la estructura terciaria** sino por el **procesamiento del ARNm** (sitio críptico de *splicing*).
5. Esto aporta una explicación estructural a Deng et al. (2025): la variante pudo alcanzar alta frecuencia **porque el cambio de aminoácido es estructuralmente tolerable**,
   dejando que la selección natural actúe sobre un defecto que está fuera de la estructura 3D.
""")
with st.expander("Proyección a futuro (mejoras posibles)"):
    st.markdown("- **Potencial electrostático superficial (APBS/Poisson–Boltzmann):** visualizar el \"parche positivo\" que crea la Lys donde antes había carga negativa.\n"
                "- **Distancias e interacciones con el hemo (PyMOL/VMD):** verificar si la Lys atrae o repele los propionatos del hemo.\n"
                "- **Explicar mecanísticamente la presión de la malaria y la talasemia** a partir del impacto de la mutación.")

# ---------------------------------------------------------------- 11
st.header("Glosario y bibliografía", anchor="s11")
with st.expander("Glosario para no especialistas"):
    G_ = {"Missense": "Mutación que cambia un aminoácido.", "Frameshift": "Inserción/deleción que corre el marco de lectura.",
          "Nonsense": "Mutación que crea un codón de parada prematuro.", "Splicing": "Corte y empalme del ARN para quitar los intrones.",
          "Promotor": "Región que regula cuánto se expresa un gen.", "Hemo": "Grupo con hierro que une el oxígeno.",
          "Estado T / R": "Formas \"tensa\" (desoxigenada) y \"relajada\" (oxigenada) de la hemoglobina.",
          "Ortólogo": "Gen equivalente en otra especie.", "Alineamiento múltiple (MSA)": "Secuencias acomodadas para comparar posición a posición.",
          "SASA": "Superficie accesible al solvente: cuánto del residuo toca el agua.", "RMSD": "Diferencia promedio de posición entre dos estructuras (en Å).",
          "Å (ángstrom)": "0,1 nanómetros.", "Factor B": "Incertidumbre posicional de un átomo en el cristal.",
          "Puente salino": "Atracción entre cargas opuestas de dos residuos cercanos.", "β⁺-talasemia": "Reducción (no ausencia) de β-globina.",
          "HbS": "Hemoglobina falciforme (β6 Glu→Val)."}
    for k, v in G_.items():
        st.markdown(f"**{k}:** {v}")
with st.expander("Bibliografía"):
    st.markdown("""
- Zhang Q, Li J, …, Xu X, Deng L. *Multi-centric origins and gene flow shape the diversity of β-thalassemia mutations in Southern East Asia.* Nat Commun. 2025;16:10220.
- Orkin SH et al. *Abnormal RNA processing due to the exon mutation of βE-globin gene.* Nature. 1982;300:768–769.
- Sen U et al. *Crystal structures of HbA2 and HbE…* Biochemistry. 2004;43:12477–12488.
- Dasgupta J et al. *Crystallization and preliminary X-ray structural studies of hemoglobin A2 and hemoglobin E…* BBRC. 2003;303:619–623.
- Salvatore F, Brunello FG, Schuster CD, Martí MA. *A structure and function-based complete mutational map of Human Hemoglobin using AI.* bioRxiv. 2025.
- Thein SL. *The molecular basis of β-thalassemia.* Cold Spring Harb Perspect Med. 2013;3(5):a011700.
- **Sección 3 (malaria):** Flatz G et al. *The 'hot-spot' of Hb E in Southeast Asia…* Hemoglobin. 2004;28:197-204 · Nagel RL et al. *Impairment of the growth of P. falciparum in HbEE erythrocytes.* J Clin Invest. 1981;68:303-305 ·
  Chotivanich K et al. (2002), vía OMIM · *Extended linkage disequilibrium surrounding the hemoglobin E variant due to malarial selection* (PMC1182083) ·
  Malaria Atlas Project, *Global Pf Parasite Rate* (PfPR2-10), versión 202608 (GBD2025), años 2000 y 2025, recortada y remuestreada a ~0,17° · Manzoni G et al. *Progress towards malaria elimination in the Greater Mekong Subregion.* Malar J. 2024 · WHO, Día Mundial del Paludismo 2021.
- **Mapa (3.5):** IthaMaps/ITHANET (prevalencia de portadores HbE por país) · Huisman et al., *A Syllabus of Human Hemoglobin Variants* (1996), ficha Hb New York β113 Val→Glu · *Clinical and molecular characteristics of hemoglobin New York in Guangxi populations*, Chin J Hematol 2013;34:696-699 · estudios de variantes de hemoglobina en Laibin (Guangxi) y norte de Guangdong (PMC10545308).
""")
    st.markdown("[UniProt P68871](https://www.uniprot.org/uniprotkb/P68871/entry) · [2HHB](https://www.rcsb.org/structure/2HHB) · "
                "[1HHO](https://www.rcsb.org/structure/1HHO) · [3DUT (HbE)](https://www.rcsb.org/structure/3DUT)")
st.divider()
st.caption("Bioinformática 2026 · Integrantes: Francisco Jose Rizzi, Nicolas Wald, Victor Orsingher.")
