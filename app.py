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


_MALARIA_B64 = (
    "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"
)  # PfPR2-10 (Malaria Atlas Project, GBD2025 v202608), 2000 y 2025, grilla de ~0,17 grados sobre lon 84-124 / lat 5-32


@st.cache_data
def cargar_malaria():
    """Capa de malaria incorporada en el propio app.py (no necesita archivos externos)."""
    import base64
    import zlib
    import numpy as np
    q = np.frombuffer(zlib.decompress(base64.b64decode("".join(_MALARIA_B64))), dtype="<u2").reshape(2, 162, 240).astype(float)
    q[q == 65535] = np.nan
    ii, jj = np.meshgrid(np.arange(162), np.arange(240), indexing="ij")
    df = pd.DataFrame({"lat": (32 - (ii * 4 + 2) * 27 / 648).ravel().round(3), "lon": (84 + (jj * 4 + 2) * 40 / 960).ravel().round(3),
                       "pfpr_2000": (q[0] / 100).ravel(), "pfpr_2025": (q[1] / 100).ravel()})
    return df.dropna().reset_index(drop=True)


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